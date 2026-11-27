"""Enumerate the lowest-layer survivors by the row whose assigned product reaches its sum."""

from argparse import ArgumentParser
from functools import lru_cache
from itertools import product
from math import gcd, isqrt, lcm, prod

from check_primitive_row_scales import direct_layers, powers, profile, scaled


@lru_cache(maxsize=None)
def divisors(n):
    factors = []
    for p, w in powers(n):
        entries = [1]
        while entries[-1] < w:
            entries.append(p*entries[-1])
        factors.append(entries)
    return sorted(prod(xs) for xs in product(*factors))


def heavy_rows(limit):
    for d in range(2, isqrt(limit//3)+1):
        for e in range(d+1, isqrt(limit+d*d)-d+1):
            if gcd(d, e) != 1:
                continue
            for f in range(e+1, (limit-d*e)//(d+e)+1):
                if gcd(d*e, f) != 1:
                    continue
                S = d*e*f
                cap = min(S, limit)
                a, b, c = e*f, d*f, d*e
                for x in range(a, cap-b-c+1, a):
                    for y in range(b, cap-x-c+1, b):
                        for z in range(c, cap-x-y+1, c):
                            if gcd(x, y, z) == 1:
                                yield tuple(sorted((x, y, z))), S


@lru_cache(maxsize=None)
def companion(v, T):
    if T == 1:
        residues = {(0, 0, 0)}
    else:
        residues = {(0, 0, 0)}
        for _, w in powers(T):
            unit = (T//w)*((v*pow(T//w, -1, w)) % w)
            residues = {tuple((x+unit*(j == i)) % T for j, x in enumerate(row))
                        for row in residues for i in range(3)}
    rows = set()
    for a, b, c in residues:
        for x in range(a or T, v-1, T):
            for y in range(b or T, v-x, T):
                z = v-x-y
                if z % T == c and gcd(x, y, z) == 1:
                    rows.add(tuple(sorted((x, y, z))))
    return rows


def enumerate_pairs(limit):
    pairs = set()
    heavy = 0
    for x, S in heavy_rows(limit):
        u = sum(x)
        heavy += 1
        for T in range((-pow(S, -1, u)) % u, isqrt(limit**3//27)+1, u):
            if gcd(S, T) != 1:
                continue
            n = S*T+1
            for v in divisors(n//u):
                if v < 3 or v > limit or gcd(u, v) != 1 or 27*T*T > v**3:
                    continue
                for y in companion(v, T):
                    first, second = sorted((x, y))
                    pairs.add((first, second, n))
    return heavy, pairs


def compare_direct(limit):
    # List directions by their entries. Grouping equal sums and H-values allows
    # the direct divisor test to compare all pairs of profiles simultaneously.
    groups = {}
    for n in range(3, limit+1):
        for a in range(1, n//3+1):
            for b in range(a, (n-a)//2+1):
                c = n-a-b
                if gcd(a, b, c) == 1:
                    H = gcd(a*b, b*c, c*a)
                    groups.setdefault((n, H), []).append((a, b, c))
    expected = set()
    for (u, h), rows in groups.items():
        if h < u:
            continue
        for (v, j), companions in groups.items():
            if h*j < u*v-1 or gcd(u, v) != 1:
                continue
            for m in divisors(lcm(h, j)):
                if (m+1) % (u*v) != 0:
                    continue
                for x, y in product(rows, companions):
                    first, second = sorted((x, y))
                    expected.add((first, second, m+1))
    assert enumerate_pairs(limit)[1] == expected
    count = sum(map(len, groups.values()))
    return count*(count+1)//2, len(expected)


def valuation(n, row, p):
    total = 0
    while n:
        n //= p
        row = tuple(a//p for a in row)
        total += n-sum(row)
    return total


def run(limit):
    heavy, pairs = enumerate_pairs(limit)
    weights = set()
    least_gcd = fewest_primes = None
    for x, y, n in pairs:
        assert profile(x, y)[2] == n
        rows = tuple(scaled(z, n) for z in (x, y))
        assert direct_layers(n, rows)
        common = []
        for p, w in powers(n)+powers(n-1):
            exponent = min(valuation(n, row, p) for row in rows)
            if exponent:
                common.append(min(w, p**exponent))
        if not common:
            raise RuntimeError((n, x, y))
        record = (prod(common), n, x, y)
        if least_gcd is None or record < least_gcd:
            least_gcd = record
        by_primes = (len(common), *record)
        if fewest_primes is None or by_primes < fewest_primes:
            fewest_primes = by_primes
        weights.add(n)
    print("Minimum endpoint gcd and witness:", least_gcd)
    print("Fewest common endpoint primes and witness:", fewest_primes)
    return heavy, len(pairs), len(weights), max(weights, default=0)


def exact_records():
    from flint import fmpz

    for n, x, y, expected in (
        (1550578, (16, 1571, 1571), (21, 141, 329), 42),
        (8055653, (44, 892, 2453), (735, 821, 821), 44),
    ):
        result = fmpz(n*(n-1))
        for row in (x, y):
            a, b, _ = scaled(row, n)
            result = result.gcd(fmpz.bin_uiui(n, a)*fmpz.bin_uiui(n-a, b))
        assert result == expected
        print("Exact witness:", n, result)


if __name__ == "__main__":
    p = ArgumentParser(description=__doc__)
    p.add_argument("--limit", type=int, default=3000)
    p.add_argument("--exact-records", action="store_true")
    args = p.parse_args()
    if args.exact_records:
        exact_records()
    print("Pairs in direct comparison:", compare_direct(200), flush=True)
    print("Heavy rows with assigned modulus, pairs, weights, largest weight:", run(args.limit))
