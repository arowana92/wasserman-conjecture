"""Check the digit estimates and count surviving primitive pairs at their sole possible weight."""

from argparse import ArgumentParser
from collections import Counter
from itertools import product
from math import comb, gcd, isqrt

from check_direction_pair_count import divisors, enumerate_pairs, heavy_rows
from check_primitive_row_scales import powers, scaled


def no_carry(p, parts):
    if p == 2:
        mask = 0
        for part in parts:
            if part & mask:
                return False
            mask |= part
        return True
    while any(parts):
        if sum(x % p for x in parts) >= p:
            return False
        parts = tuple(x // p for x in parts)
    return True


def order_at(n, p):
    value = 0
    while n % p == 0:
        value += 1
        n //= p
    return value


def chain_count(p, length):
    # Start with a fixed zero before the chain and allow a+b < p at each edge.
    current = [1] * p
    for _ in range(length):
        current = [sum(current[:p-i]) for i in range(p)]
    return current[0]


def digit_controls():
    trials = tight = balanced = 0
    for p in (2, 3, 5, 7):
        for width in range(1, 5):
            size = p**width
            chains = [chain_count(p, j) for j in range(width+1)]
            for gap in range(1, width+2):
                bound = 1
                for i in range(min(gap, width)):
                    bound *= chains[1+(width-1-i)//gap]
                # An arbitrary step, low translated digits, and different unit slopes.
                for start, step, a, b in (
                        (0, 1, 1, 1), (1, p, 1, 1), (37, 2*p, p+1, 2*p+1),
                        (13*p+1, p**2+1, p**2-1, p**2+1)):
                    count = sum(no_carry(p, (a*w, b*p**gap*w))
                                for w in range(start, start+step*size, step))
                    assert count <= bound
                    if start == 0:
                        assert count == bound
                        tight += 1
                    trials += 1
            for a, b in product(range(1, p), repeat=2):
                # The current two output digits range through permutations.
                for c, d in product(range(p), repeat=2):
                    count = sum((a*j+c) % p + (b*j+d) % p < p for j in range(p))
                    if (a+b) % p:
                        assert 2*count == p+1
                        balanced += 1
    print(f"Digit checks: {trials} progressions, {tight} equalities, "
          f"{balanced} equal-delay columns.")


def heavy_crosscheck(limit, cap):
    # Reconstruct the heavy row from two congruences for W = multiplier*f.
    # Compare with traversing the heavy primitive entries, then trying multipliers.
    reconstructed = set()
    for d in range(2, isqrt(limit//3)+1):
        for e in range(d+1, isqrt(limit+d*d)-d+1):
            if gcd(d, e) != 1:
                continue
            de = d*e
            for a in range(1, d):
                for b in range(1, e):
                    slopes = a*e, b*d
                    total = sum(slopes)
                    bound_f = (limit-de)//total
                    if bound_f <= e or gcd(slopes[0], d)*gcd(slopes[1], e) != 1:
                        continue
                    start = pow(slopes[0], -1, d)
                    start += d*((pow(slopes[1], -1, e)-start)*pow(d, -1, e) % e)
                    for w in range(start, cap//total+1, de):
                        for f in divisors(w):
                            if not e < f <= bound_f or gcd(de, f) > 1:
                                continue
                            m = w//f
                            if gcd(m, f) > 1:
                                continue
                            z = pow(m*de, -1, f)
                            row = tuple(sorted((a*e*f, b*d*f, z*de)))
                            u, s = sum(row), de*f
                            n = m*u
                            if (u <= min(limit, s) and n <= cap and gcd(*row) == 1
                                    and (n-1) % s == 0 and gcd(s, (n-1)//s) == 1):
                                reconstructed.add((row, s, m))
    direct = set()
    for row, s in heavy_rows(limit):
        u = sum(row)
        for m in range(pow(u, -1, s), cap//u+1, s):
            if gcd(s, (m*u-1)//s) == 1:
                direct.add((row, s, m))
    assert direct == reconstructed and direct
    print(f"Both reconstructions give {len(direct):,} heavy-row assignments.")


def progression_rows():
    # Compare the light-row block count with every positive row in each CRT class.
    checks = 0
    for p, v, step, m in product((2, 3, 5), range(3, 44), range(1, 14), (1, 2, 3, 6, 15)):
        counts = Counter()
        for a in range(1, v-1):
            for b in range(1, v-a):
                row = (a, b, v-a-b)
                if no_carry(p, tuple(m*x for x in row)):
                    counts[tuple(x % step for x in row)] += 1
        # Strictly more than the greatest possible quotient difference.
        width, span = 0, 1
        while span <= (v-3)//step:
            width += 1
            span *= p
        n = (m*v)//p**order_at(m*step, p)
        bound = 1
        for _ in range(width):
            bound *= comb(n % p+2, 2)
            n //= p
        assert all(count <= bound for count in counts.values())
        checks += len(counts)
    print(f"Light-row digit bound checked on {checks:,} nonempty residue classes.")


def census(limit):
    _, pairs = enumerate_pairs(limit)
    prime_list = (2, 3, 5, 7)
    avoidable = Counter()
    survivors = []
    common_sets = set()
    for x, y, n in pairs:
        rows = (scaled(x, n), scaled(y, n))
        avoids = {p for p in prime_list if any(no_carry(p, row) for row in rows)}
        avoidable.update(avoids)
        if 2 in avoids:
            common = tuple(sorted(p for p, _ in powers(n)+powers(n-1)
                                  if all(not no_carry(p, row) for row in rows)))
            assert common
            common_sets.add(common)
            survivors.append((max(sum(x), sum(y)), n, x, y, common))
    # A fixed endpoint prime together with 2 can certify the entire finite list
    # precisely if it lies in every list attached to an odd-coefficient pair.
    shared = sorted(set.intersection(*(set(s) for s in common_sets))) if common_sets else []
    print(f"Through {limit}, {len(pairs):,} lowest-layer pairs. Counts with an "
          f"avoiding row: {[(p, avoidable[p]) for p in prime_list]}.")
    print("Common endpoint primes over the parity survivors:", shared)
    if survivors:
        print("First by primitive sum:", min(survivors))
    if limit >= 2769:
        assert min(survivors)[:4] == (2769, 3126201, (49, 296, 784), (275, 319, 2175))
    if limit == 3000:
        assert len(pairs) == 75278 and [avoidable[p] for p in prime_list] == [1, 7, 12, 31]
    if limit == 6000:
        assert len(pairs) == 551062 and [avoidable[p] for p in prime_list] == [5, 11, 41, 109]
        assert shared == [3, 5]
    if limit >= 400:
        for first, second, n in (
                ((5, 31, 31), (33, 78, 286), 26599),
                ((9, 127, 127), (95, 120, 152), 96521)):
            assert (first, second, n) in pairs
        assert gcd(26599*26598, 96521*96520) == 2


def exact_witness():
    from flint import fmpz

    n = 3126201
    value = fmpz(n*(n-1))
    for row in ((49, 296, 784), (275, 319, 2175)):
        a, b, _ = scaled(row, n)
        k = fmpz.bin_uiui(n, a)*fmpz.bin_uiui(n-a, b)
        assert bool(k % 2) == no_carry(2, scaled(row, n))
        value = value.gcd(k)
    assert value == 1221641195775
    print("Exact gcd with the endpoint at the first parity survivor:", value)


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--limit", type=int, default=3000)
    parser.add_argument("--exact", action="store_true", help="use python-flint at the first witness")
    args = parser.parse_args()
    digit_controls()
    heavy_crosscheck(140, 30000)
    progression_rows()
    census(args.limit)
    if args.exact:
        exact_witness()
