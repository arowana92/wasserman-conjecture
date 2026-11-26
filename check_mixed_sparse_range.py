"""Exact lower range for one binary and one odd sparse predecessor prime."""

from argparse import ArgumentParser
from bisect import bisect_left, bisect_right
from collections import Counter
from itertools import combinations
from math import gcd, isqrt

from check_dense_offsets import free_rows


def next_sparse(n):
    """Smallest integer at least n containing at most six ones."""
    tail = n
    for _ in range(6):
        if not tail:
            return n
        bit = 1 << (tail.bit_length()-1)
        tail -= bit
    return n-tail+bit if tail else n


def search(p, limit):
    powers = [p]
    while powers[-1]*p <= limit:
        powers.append(powers[-1]*p)
    hits = set()

    def visit(total, last, left, copies):
        target = next_sparse(total+left*p)
        if target > limit or target > total+left*powers[last]:
            return
        lo = bisect_left(powers, (target-total+left-1)//left, 0, last+1)
        hi = bisect_right(powers, limit-total-(left-1)*p, lo, last+1)
        for j in range(lo, hi):
            used = copies+1 if j == last else 1
            if used >= p:
                continue
            n = total+powers[j]
            if left > 1:
                visit(n, j, left-1, used)
            elif n == next_sparse(n) and bin(n).count('1') in (5, 6):
                hits.add((n+1, p))

    visit(0, len(powers)-1, 6, 0)
    return hits


def binary_cross_check(primes, bits):
    hits = set()
    for count in (5, 6):
        for positions in combinations(range(1, bits), count):
            n = sum(1 << e for e in positions)
            for p in primes:
                if n % p:
                    continue
                total = 0
                remaining = n
                while remaining and total <= 6:
                    remaining, digit = divmod(remaining, p)
                    total += digit
                if total == 6:
                    hits.add((n+1, p))
    return hits


def check(bits):
    primes = [p for p in range(3, 256, 2) if all(p%r for r in range(2, isqrt(p)+1))]
    hits = set()
    for p in primes:
        hits.update(search(p, 2**bits-2))
    assert hits == binary_cross_check(primes, 34)
    counts = Counter()
    for n, p in hits:
        m = n-1
        while m % p == 0:
            m //= p
        for row in free_rows(n, p):
            counts['rows'] += 1
            u = gcd(*row)
            if u < 7:
                continue
            counts['large row gcds'] += 1
            s = 1
            for x, y in combinations(row, 2):
                s *= gcd(m, x, y)
            s //= s & -s
            assert 27*m*m > u**3*s*s, (n, p, row)
    assert len(hits) == 1987 and max(hits) == (10737942553, 3)
    print(f'{len(hits)} weight/base pairs, largest {max(hits)}, {dict(counts)}.')


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--bits', type=int, default=128, help='exact ceiling exponent, at least 34')
    args = parser.parse_args()
    if args.bits < 34:
        parser.error('--bits must be at least 34 for the independent binary comparison')
    check(args.bits)
