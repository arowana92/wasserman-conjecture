"""Check the gcd inequality by exact coefficients and decide the finite survivors."""

from argparse import ArgumentParser
from collections import Counter
from functools import lru_cache
from itertools import product
from math import comb, gcd, isqrt
from pathlib import Path

from check_endpoint_pair_bound import allocate, prime_powers


def powers(p, bound):
    value = 1
    while value <= bound:
        yield value
        value *= p


@lru_cache(maxsize=None)
def row_pairs(n, p):
    # Symmetric rows suffice. Avoidance of p is tested from the coefficients.
    values = [comb(n, a) * comb(n - a, b) for a in range(1, n // 3 + 1)
              for b in range(a, (n - a) // 2 + 1)]
    weight = tuple(prime_powers(n))
    predecessor = tuple(prime_powers(n - 1))
    distribution = Counter()
    for first in values:
        if first % p == 0:
            continue
        for second in values:
            u, _, partial0 = allocate(weight, first, second)
            s, _, partial1 = allocate(predecessor, first, second)
            g = gcd(n * (n - 1), first, second)
            distribution[u**3 * s**2, g, partial0 or partial1] += 1
    return distribution


def exact_check(limit=48):
    counts = Counter()
    for p in range(2, limit + 1):
        if any(p % k == 0 for k in range(2, isqrt(p) + 1)):
            continue
        shifts = list(powers(p, limit))
        for high, q, low in product(shifts[1:], shifts[1:], shifts[1:]):
            if high < low:
                continue
            for m in range(1, (limit - high - low) // q + 1):
                n = high + m * q + low
                bound = (low * m)**5
                for (assigned, g, partial), count in row_pairs(n, p).items():
                    assert assigned <= bound
                    assert 27 * (n - 1)**2 <= g**3 * bound
                    counts["coefficient pairs"] += count
                    if partial:
                        counts["partially retained"] += count
                    if 27 * (n - 1)**2 > bound:
                        counts["forcing a divisor"] += count
                    if high == low and q > high and p == 2:
                        counts["coincident binary powers"] += count
    return counts


def cross_check_rows(path, cap, min_low):
    # Independent candidate generation on a subset of the coefficient range.
    # Use ordinary ordered digits and test every k from 1 through m.
    mode = path.stem.split("two-singleton-")[-1]
    if mode not in {"binary-ordered", "binary-above", "odd"}:
        raise ValueError("expected two-singleton-{binary-ordered,binary-above,odd}.rows")
    recorded = {tuple(map(int, row.split())) for row in path.read_text().splitlines()}
    recorded = {r for r in recorded if r[1]*r[4] <= cap and r[4] >= min_low}
    computed = set()
    weights = 0
    def fits(x, y, prime):
        while x:
            if x % prime > y % prime:
                return False
            x //= prime
            y //= prime
        return True
    for p in range(2, cap + 1):
        if (p == 2) == (mode == "odd") or any(p % z == 0 for z in range(2, isqrt(p)+1)):
            continue
        for low in list(powers(p, cap))[1:]:
            if low < min_low:
                continue
            for m in range(1, cap // low + 1):
                top = isqrt((low*m)**5 // 27) + 1
                if top <= 10**9:
                    continue
                for high in powers(p, top):
                    if high < low:
                        continue
                    for q in powers(p, (top - high - low) // m):
                        if mode == "binary-above" and q <= high or mode == "binary-ordered" and q > high:
                            continue
                        n = high + m*q + low
                        if n <= 10**9 or n*(n-1) % p:
                            continue
                        weights += 1
                        target = 27*(n-1)**2
                        if q > 2*low and high == low and p == 2 and low**3*m**4 < 2*target:
                            continue
                        # The fifth-power consequence also gives a cheap threshold
                        # for the brute-force list of individual parts.
                        if q <= low and (m*q)**5 <= target:
                            continue
                        for k in range(1, m+1):
                            z = k*q
                            if z < 3 or n-z < 6 or not fits(k,m,p) or not fits(z,n,p):
                                continue
                            u = gcd(n,z)
                            s = gcd(n-1,z*(z-1))
                            if u < 3 or u**3*s**2 < target:
                                continue
                            if q <= low and z**5 <= target:
                                continue
                            if q > low and (high > low or p != 2):
                                if low**3*gcd(n//low,k)*(k*(m+(low-1)*k))**2 < target:
                                    continue
                            computed.add((p,m,high,q,low,n))
                            break
    assert computed == recorded, (computed-recorded, recorded-computed)
    return weights, len(computed)


def finite(paths):
    from check_binary_gap_finite import check_weight
    counts = Counter()
    seen = set()
    for path in paths:
        rows = [tuple(map(int, line.split())) for line in path.read_text().splitlines()]
        for p, m, high, q, low, n in rows:
            assert n == high + m*q + low and n*(n-1) % p == 0
            assert all(v in powers(p, v) for v in (high, q, low))
        seen.update(row[-1] for row in rows)
        print(path.name, len(rows), "rows,", len({row[-1] for row in rows}), "weights")
    for n in sorted(seen):
        check_weight(n, counts)
    print(len(seen), "distinct weights,", dict(counts))


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("rows", nargs="*", type=Path)
    parser.add_argument("--cross-cap", type=int, default=0)
    parser.add_argument("--cross-min-low", type=int, default=2)
    args = parser.parse_args()
    print("Exact checks:", dict(exact_check()))
    if args.cross_cap:
        for path in args.rows:
            print("Independent scan:", path.name, cross_check_rows(path, args.cross_cap, args.cross_min_low))
    if args.rows:
        finite(args.rows)
