"""Reverse the cross-part search via the actual (unnormalized) parts."""

from argparse import ArgumentParser
from collections import Counter
from functools import lru_cache
from itertools import product
from math import comb, gcd, prod
import re

from search_even_endpoint import avoids, prime_powers


@lru_cache(maxsize=100_000)
def factors(n):
    return prime_powers(n)


def divisors(n):
    values = [1]
    for p, power in factors(n).items():
        old = values.copy()
        while power > 1:
            values.extend(d * power for d in old)
            power //= p
    return values


def moment_pairs(t):
    for k in range(1, t):
        target = t * (t - k) * (t*t - k)
        for first in divisors(target):
            second = target // first
            if (first + t*t) % k or (second + t*t) % k:
                continue
            a, b = (first + t*t) // k, (second + t*t) // k
            if a*b % t == 0:
                yield a, b, a*b // t


def assignments(n, a, b):
    # Assign each full predecessor power to one row using actual residues.
    choices = []
    for w in factors(n-1).values():
        positions = [i for i, value in enumerate((a, a-1, b, b-1)) if value % w == 0]
        if not positions:
            return
        choices.append([(i, w) for i in positions])
    for assignment in product(*choices):
        yield tuple(prod(w for i, w in assignment if i == j) for j in range(4))


def normalized_splits(distinguished, modulus, multiplier, total):
    spare = total - distinguished
    if spare <= 0 or spare % multiplier:
        return []
    s = spare // multiplier
    if gcd(s, modulus) != 1:
        return []
    return [(distinguished, i*multiplier, (s-i)*multiplier)
            for i in range(1, s//2+1) if i*(s-i) % modulus == 0]


def check_moments(limit, first_t=1):
    counts = Counter()
    largest = 0
    for t in range(first_t, limit + 1):
        for a, b, n in moment_pairs(t):
            for x in divisors(t):
                y = t // x
                if a % x or b % y:
                    continue
                u, v = a // x, b // y
                if min(u, v) < 3 or gcd(u, v) > 1:
                    continue
                for f, left, h, right in assignments(n, a, b):
                    if f >= h or x % f or y % h:
                        continue
                    rows_a = normalized_splits(x, f, left, v)
                    rows_c = normalized_splits(y, h, right, u)
                    if not rows_a or not rows_c:
                        continue
                    counts['systems'] += 1
                    counts['pairs'] += len(rows_a) * len(rows_c)
                    largest = max(largest, n)
                    for row in rows_c:
                        actual = tuple(v*i for i in row)
                        if all(avoids(n, actual, p) for p in factors(v*h*right)):
                            counts['passing second rows'] += 1
                            for row_a in rows_a:
                                actual_a = tuple(u*i for i in row_a)
                                if all(avoids(n, actual_a, p) for p in factors(u*f*left)):
                                    raise ArithmeticError((n, actual_a, actual))
    interval = f"Through {limit}" if first_t == 1 else f"From {first_t} through {limit}"
    print(f"{interval}: {counts['systems']} systems, {counts['pairs']} row pairs, "
          f"{counts['passing second rows']} passing second rows, 0 passing pairs. "
          f"Largest weight {largest}")


def audit_pairs(path):
    checked = 0
    weights = set()
    second_rows = set()
    expected = None
    largest_witness = 0
    least_binary_valuation = None
    with open(path) as records:
        for line in records:
            if line.startswith('Through '):
                summary = re.search(r'(\d+) passing second rows, (\d+) passing pairs', line)
                assert summary and int(summary[2]) == 0 and expected is None
                expected = int(summary[1])
                continue
            if not line.startswith('Pair '):
                raise ValueError(line)
            assert expected is None
            t, n, u, s, v, h, *parts = map(int, line.split()[1:])
            a, c = tuple(parts[:3]), tuple(parts[3:])
            assert len(c) == 3 and min(parts) > 0 and sum(a) == sum(c) == n
            assert n == u*v and n-1 == s*h and n*t == a[0]*c[0]
            assert all(avoids(n, c, p) for p in factors(v*h))
            common = [p for p in (*factors(n), *factors(n-1))
                      if not avoids(n, c, p) and not avoids(n, a, p)]
            assert common, (n, a, c)
            largest_witness = max(largest_witness, min(common))
            binary_valuation = min(sum(bin(x).count('1') for x in row)
                                   - bin(n).count('1') for row in (a, c))
            if least_binary_valuation is None or binary_valuation < least_binary_valuation:
                least_binary_valuation = binary_valuation
            checked += 1
            weights.add(n)
            second_rows.add((n, u, s, v, h, a[0], c))
    assert expected == len(second_rows), (expected, len(second_rows))
    print(f"Audited {checked} pairs at {len(weights)} weights. "
          f"Largest least common endpoint prime: {largest_witness}. "
          f"Least common binary valuation of the two rows: {least_binary_valuation}")


def check_coefficients(limit):
    comparisons = 0
    for n in range(3, limit + 1):
        rows = [((a, b, n-a-b), comb(n, a) * comb(n-a, b))
                for a in range(1, n//3+1) for b in range(a, (n-a)//2+1)]
        for (left, left_k), (right, right_k) in product(rows, repeat=2):
            common = gcd(left_k, right_k)
            d, j = gcd(common, n), gcd(common, n-1)
            for i, a in enumerate(left):
                for k, b in enumerate(right):
                    other_four = prod(left[:i] + left[i+1:] + right[:k] + right[k+1:])
                    assert (other_four*d*d*j) % (n*n*(n-1)) == 0
                    assert a*b*d % n == 0
                    if min(a, b) >= 2:
                        assert 9*d**13*j**6*(a*b)**9 >= n**9*(n-1)**2
                    comparisons += 1
    print(f"Exact trinomial comparisons through weight {limit}: {comparisons}")


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--limit', type=int, default=200)
    parser.add_argument('--first', type=int, default=1)
    parser.add_argument('--weight', type=int, default=35)
    parser.add_argument('--audit', help='audit the surviving second rows from the C++ search')
    args = parser.parse_args()
    if args.audit:
        audit_pairs(args.audit)
    else:
        check_coefficients(args.weight)
        check_moments(args.limit, args.first)
