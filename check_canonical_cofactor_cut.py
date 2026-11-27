"""Check the inverse residues, canonical cut and descending digit chains."""

from argparse import ArgumentParser
from bisect import bisect_right
from math import comb, gcd, isqrt

from check_three_factors import powers, valuation
from prime_certificates import free_count, free_rows

from check_static_power_transfer import prime_powers
from check_two_power_seed import free, grouped, split_costs


def canonical(a, b):
    c, d = pow(b, -1, a), pow(a, -1, b)
    assert b*c + a*d == a*b + 1
    return c, d, gcd(c, d)


def necessary(a, b, p, q):
    c, d, g = canonical(a, b)
    return g >= 3 and free((0, d, b-d), p) and free((0, c, a-c), q)


def descending(u, p, i):
    digits = []
    while u:
        u, digit = divmod(u, p)
        digits.append(digit)
    digits += [0]*i
    return digits[i] < digits[0] and all(digits[j] <= digits[j-i]
                                          for j in range(i+1, len(digits)))


def ordinary_prime_example():
    n = 2_100_875_036
    w = 103 * 582_767
    assert powers(w) == {103: 103, 582767: 582767}
    assert n == 35*w + 1 and (w-1) % (5**5 * 7**4) == 0
    assert not necessary(5, 7, 5, 7)
    primes = (7, 563, 71761, n-3)
    assert all(all(p % d for d in range(2, isqrt(p)+1)) for p in primes)
    # Disjoint avoidance sets suffice for four candidate primes and three rows.
    for p, others in ((7, primes[1:]), (n-3, primes[1:3])):
        ordered = 0
        for row in free_rows(n, p):
            assert valuation(n, row, p) == 0
            assert all(valuation(n, row, q) > 0 for q in others)
            ordered += (1, 3, 6)[len(set(row))-1]
        assert ordered == free_count(n, p)
    unit = 563 * 71761
    total = n // unit
    assert n == unit*total and total == 52
    for a in range(1, total//3 + 1):
        for b in range(a, (total-a)//2 + 1):
            row = (a*unit, b*unit, (total-a-b)*unit)
            assert any(valuation(n, row, q) > 0 for q in primes[1:3])
    print(f"direct four-prime certificate at {n}: {primes}")


def companion_projection():
    checked = nonprimitive = 0
    for p in (2, 3, 5, 7, 11):
        for n in range(4, 81):
            if n % p == 0:
                continue
            layer = p
            while layer <= n*n:
                layer *= p
            for u in range(3, n):
                if u % p == 0:
                    continue
                v = n * pow(u, -1, layer) % layer
                for f in range(1, u//3+1):
                    for g in range(f, (u-f)//2+1):
                        parts = (f, g, u-f-g)
                        residues = tuple(v*x % layer for x in parts)
                        if sum(residues) != n:
                            continue
                        h = gcd(f, gcd(g, u-f-g))
                        assert n % (u//h) == 0
                        assert all(u*r == n*x for r, x in zip(residues, parts))
                        checked += 1
                        nonprimitive += (n % u != 0)
    print(f"{checked} bounded projections, including {nonprimitive} with scale not dividing seed")
    assert nonprimitive


def check(limit, cross):
    powers = prime_powers(limit)
    candidates = survivors = general = 0
    values = sorted(powers)
    for index, a in enumerate(values):
        if a*a >= limit:
            break
        p = powers[a]
        for j in range(index + 1, bisect_right(values, limit//a)):
            b = values[j]
            q = powers[b]
            if p == q:
                continue
            c, d, g = canonical(a, b)
            n = a*b + 1
            ok = necessary(a, b, p, q)
            if a*b <= cross:
                shared = comb(n, b*c)
                assert (shared % p != 0) == free((0, d, b-d), p)
                assert (shared % q != 0) == free((0, c, a-c), q)
                for x, u, _, _ in grouped(a*b):
                    assert x in (b*c, a*d) and g % u == 0 and ok
                if not ok:
                    assert split_costs(a*b)[1] >= a*b, (a, b, c, d, g)
                    general += 1
            candidates += 1
            survivors += ok
    digit_cases = 0
    for p, i in ((2, 2), (2, 3), (2, 4), (3, 1), (3, 2), (3, 3), (5, 1), (5, 2), (7, 1)):
        a = p**i
        for b in range(a-1, 100001, a):
            u = b+1
            while u % p == 0:
                u //= p
            assert descending(u, p, i) == free((0, (b+1)//a, b-(b+1)//a), p)
            digit_cases += 1
    print(f"{candidates} unordered cofactor pairs through product {limit}, "
          f"{survivors} survive the canonical test")
    print(f"{general} exclusions cross-checked through product {cross}")
    print(f"{digit_cases} descending digit tests, including non-prime-power cofactors")


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--limit", type=int, default=1_000_000)
    parser.add_argument("--cross", type=int, default=500)
    args = parser.parse_args()
    if not 1 <= args.cross <= args.limit:
        parser.error("require 1 <= cross <= limit")
    check(args.limit, args.cross)
    ordinary_prime_example()
    companion_projection()
