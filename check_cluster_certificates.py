"""Enumerate common four-prime lists for pairs d*P*Q*R and d*P*Q*R+1."""

from argparse import ArgumentParser
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from math import prod

from flint import fmpz


PRIMES = (1000000007, 1000000009, 1000000021)
PRODUCT = prod(PRIMES)
LIMIT = PRIMES[0] // (2*(PRIMES[2]-PRIMES[0])**2)


def is_prime(n):
    return bool(fmpz(n).is_prime())


def universal_prime(n):
    """Return a prime at which n has at most two digit units, or zero."""
    if bin(n).count('1') <= 2:
        return 2
    if n % 2:
        if is_prime(n):
            return n
        # In an odd base an odd digit sum of at most two is one.
        if not fmpz(n).is_perfect_power():
            return 0
        for i in range(2, n.bit_length()):
            b = int(fmpz(n).root(i))
            if b**i == n and is_prime(b):
                return b
        return 0
    # For odd b, the two units either coincide or have different positions.
    half = n//2
    if is_prime(half):
        return half
    if is_prime(n-1):
        return n-1
    if fmpz(half).is_perfect_power():
        for i in range(2, half.bit_length()):
            b = int(fmpz(half).root(i))
            if b**i == half and is_prime(b):
                return b
    value = fmpz(n)
    for i in range(2, n.bit_length()):
        b = int(value.root(i))
        if b <= 2:
            break
        rest = n-b**i
        if rest <= 0:
            continue
        while rest % b == 0:
            rest //= b
        if rest == 1 and is_prime(b):
            return b
    return 0


def certificates(bounds):
    lo, hi = bounds
    counts = Counter()
    max_gap = (-1, 0)
    for d in range(lo, hi+1):
        n = d*PRODUCT
        fourth = n if n%2 else n-1
        while not is_prime(fourth):
            fourth -= 2
        gap = n-fourth
        assert 2*fourth > n+1 and gap+1 <= min(PRIMES)
        max_gap = max(max_gap, (gap, d))
        counts[bool(universal_prime(n)), bool(universal_prime(n+1))] += 1
    return counts, max_gap


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--first', type=int, default=1)
    parser.add_argument('--limit', type=int, default=LIMIT)
    parser.add_argument('--jobs', type=int, default=1)
    args = parser.parse_args()
    if not 1 <= args.first <= args.limit <= LIMIT or args.jobs < 1:
        parser.error(f'multiplier bounds must be in [1,{LIMIT}] and jobs must be positive')
    p, q, r = PRIMES
    assert all(is_prime(v) for v in PRIMES) and r-q > q-p
    blocks = [(i, min(i+8191, args.limit)) for i in range(args.first, args.limit+1, 8192)]
    counts = Counter()
    largest = (-1, 0)
    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        for block_counts, gap in pool.map(certificates, blocks):
            counts.update(block_counts)
            largest = max(largest, gap)
    assert sum(counts.values()) == args.limit-args.first+1
    print(f'Multipliers {args.first} through {args.limit}; largest near-prime gap {largest}.')
    for a in (False, True):
        for b in (False, True):
            print(f'Minimum sizes ({1 if a else 4}, {1 if b else 4}): {counts[a,b]}')


if __name__ == '__main__':
    main()
