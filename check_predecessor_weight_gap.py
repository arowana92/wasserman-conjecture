"""Reproduce the weight-gap tests with the Python standard library."""

from argparse import ArgumentParser
from collections import Counter
from itertools import combinations_with_replacement
from math import comb, gcd, isqrt, prod


def factorizations(lo, hi):
    """Full prime powers for each integer in the closed positive interval."""
    remainders = list(range(lo, hi + 1))
    rows = [{} for _ in remainders]
    limit = isqrt(hi)
    prime = bytearray(b"\x01") * (limit + 1)
    for p in range(2, limit + 1):
        if not prime[p]:
            continue
        if p * p <= limit:
            prime[p*p::p] = bytes((limit - p*p) // p + 1)
        for i in range(-lo % p, len(rows), p):
            power = 1
            while remainders[i] % p == 0:
                remainders[i] //= p
                power *= p
            rows[i][p] = power
    for i, r in enumerate(remainders):
        if r > 1:
            rows[i][r] = r
        assert prod(rows[i].values()) == lo + i
    return rows


def criterion(n, weight, predecessor):
    if len(predecessor) < 2 or len(weight) < 2:
        return "outside", None
    q, r = sorted(predecessor.values())[-2:]
    w = q * r
    d = (n - 1) // w
    gaps = []
    forced = []
    for p, power in weight.items():
        tail = n // power - 1
        b = 1
        while tail % p == 0:
            b *= p
            tail //= p
        gaps.append(b)
        inverse = pow(w, -1, b) if b > 1 else 1
        if power * inverse >= d:
            forced.append(b)
    if max(gaps) >= d:
        return "one gap", d
    if d > 16 * prod(gaps):
        return "inconclusive", None
    if not forced:
        return "inconclusive", None
    forced.sort(reverse=True)
    if 2 * forced[0] + 1 >= d:
        return "one forced prime", d
    if len(forced) >= 2 and forced[0] * forced[1] >= d:
        return "two forced primes", d
    if 27 * prod(forced)**2 > (d-1)**3:
        return "product of forced gaps", d
    return "inconclusive", None


def exact_check(limit):
    factors = factorizations(2, limit)
    certified = row_count = pair_count = 0
    for n in range(3, limit + 1):
        tag, _ = criterion(n, factors[n-2], factors[n-3])
        if tag in ("outside", "inconclusive"):
            continue
        certified += 1
        endpoint = n * (n - 1)
        divisors = set()
        for a in range(1, n // 3 + 1):
            choose_a = comb(n, a)
            for b in range(a, (n-a)//2 + 1):
                k = choose_a * comb(n-a, b)
                divisors.add(gcd(endpoint, k))
                row_count += 1
        for x, y in combinations_with_replacement(divisors, 2):
            assert gcd(x, y) > 1, (n, x, y)
            pair_count += 1
    print(f"Exact coefficients through {limit}: {certified} certified weights, "
          f"{row_count} rows, {pair_count} endpoint-gcd pairs")


def check_progressions():
    tests = 0
    for m in range(1, 8):
        for r in (5, 13, 29, 37):
            for c in range(1, 64, 2):
                v = 1 if c % 4 == 1 else 2
                modulus = 2**(m+v)
                period = 2**(m+v-2)
                residues = [b for b in range(period)
                            if (1+c*pow(r, b, modulus)) % modulus == 2**v]
                assert len(residues) == 1 and period <= 2**m
                tests += 1
    print(f"Binary progression checks: {tests}")


def census(start, count):
    factors = factorizations(start-1, start+count-1)
    totals = Counter()
    largest = {}
    for i in range(count):
        n = start + i
        tag, d = criterion(n, factors[i+1], factors[i])
        totals[tag] += 1
        if d is not None and (tag not in largest or d > largest[tag][0]):
            largest[tag] = (d, n)
    print(f"Census [{start}, {start+count-1}]: {dict(sorted(totals.items()))}")
    print(f"Largest cofactor and its first weight by certificate: {largest}")


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--exact", type=int, default=300)
    parser.add_argument("--start", type=int, default=10**12 + 1234567890)
    parser.add_argument("--count", type=int, default=100000)
    args = parser.parse_args()
    if args.exact < 3 or args.start < 3 or args.count < 1:
        parser.error("need exact and start >= 3 and count >= 1")
    exact_check(args.exact)
    check_progressions()
    census(args.start, args.count)
