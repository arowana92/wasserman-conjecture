"""Independent finite arithmetic checks for six units at a weight prime."""

from argparse import ArgumentParser
from collections import Counter
from functools import lru_cache
from itertools import combinations
from math import comb, gcd, prod

from flint import fmpz
from prime_certificates import free_count, free_in_base
from search_endpoint import smallest_factors


def direct_coefficients(limit):
    spf = smallest_factors(limit)
    cases = 0
    least = None
    for n in range(3, limit + 1):
        bases = []
        for p in range(2, n + 1):
            if spf[p] != p or n % p or not free_count(n, p):
                continue
            value, digit_sum = n, 0
            while value:
                value, digit = divmod(value, p)
                digit_sum += digit
            if digit_sum <= 6:
                bases.append(p)
        for a in range(1, n // 3 + 1):
            for b in range(a, (n - a) // 2 + 1):
                coefficient = comb(n, a) * comb(n - a, b)
                for p in bases:
                    if coefficient % p:
                        endpoint_gcd = gcd(coefficient, n * (n - 1))
                        assert endpoint_gcd > 2 * n
                        cases += 1
                        if least is None or endpoint_gcd * least[1] < least[0] * n:
                            least = endpoint_gcd, n, (a, b, n - a - b), p
    print('Exact sparse-row gcd comparisons', cases, 'smallest observed ratio witness', least)


def brute(bits):
    counts = Counter()
    for extras in combinations(range(2, bits), 5):
        units = (2,) + tuple(1 << i for i in extras)
        n = sum(units)
        m = n - 1
        counts['weights'] += 1
        for j in range(1, 6):
            rest = [x for i, x in enumerate(units) if i not in (0, j)]
            for k in range(1, 4):
                pairs = [(units[0], units[j]), (rest[0], rest[k]),
                         tuple(x for i, x in enumerate(rest) if i not in (0, k))]
                a, b, c = map(sum, pairs)
                u = gcd(a, b, c)
                s = gcd(m, a - 1, b, c) * gcd(m, a, b - 1, c) * gcd(m, a, b, c - 1)
                counts['pairings'] += 1
                if u >= 3 and 27 * m * m <= u**3 * s**2:
                    counts['passed size'] += 1
                    assert s == 1 and len({y // x for x, y in pairs}) == 1
    print('Direct six-bit enumeration', bits, dict(counts))


@lru_cache(maxsize=None)
def factored(n):
    factors = {int(p): int(p**e) for p, e in fmpz(n).factor()}
    assert prod(factors.values()) == n
    assert all(fmpz(p).is_prime() for p in factors)
    return factors


def independent_blocks(bits):
    counts = Counter()
    for middle, high in combinations(range(2, bits - 3), 2):
        head = 2 + (1 << middle) + (1 << high)
        for shift in range(high + 1, bits - high):
            row = tuple(x * ((1 << shift) + 1) for x in (2, 1 << middle, 1 << high))
            n = sum(row)
            m = n - 1
            if 27 * m * m > gcd(*row)**3:
                continue
            counts['size candidates'] += 1
            # Test every prime in the entire weight to find the mandatory scale.
            factors = factored(head).copy()
            for p, w in factored((1 << shift) + 1).items():
                factors[p] = factors.get(p, 1) * w
            assert prod(factors.values()) == n
            scale = prod(w for p, w in factors.items() if not free_in_base(row, p))
            u = n // scale
            if 27 * m * m > u**3:
                continue
            counts['binary row passes size'] += 1
            atoms = factored(m)
            if len(atoms) < 3 or any(2*w+1 > u for w in atoms.values()):
                continue
            counts['companion CRT systems'] += 1
            groups = [1]
            for w in atoms.values():
                groups += [w * v for v in groups]
            entries = {}
            for group in groups:
                if group in (1, m):
                    continue
                other = m // group
                multiplier = pow(other, -1, group)
                if multiplier % scale == 0:
                    entries[group] = other * multiplier
            counts['eligible entries'] += len(entries)
            for a, b, c in combinations(entries, 3):
                if a * b * c == m and gcd(a, b) == gcd(a*b, c) == 1:
                    actual = sorted(entries[g] for g in (a, b, c))
                    if sum(actual) == n:
                        counts['companion triples'] += 1
                        assert not all(free_in_base(actual, p) for p in {*atoms, *factored(scale)})
    print('Independent predecessor partitions', bits, dict(counts))


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--bits', type=int, default=24)
    parser.add_argument('--blocks', type=int, default=160)
    args = parser.parse_args()
    direct_coefficients(160)
    brute(args.bits)
    independent_blocks(args.blocks)
