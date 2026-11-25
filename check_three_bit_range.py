"""Certify all three-row selections at N with exactly three binary ones below 2**bits."""

from argparse import ArgumentParser
from collections import Counter
from itertools import combinations
from math import gcd, lcm

from flint import fmpz

from search_even_endpoint import avoids, prime_powers


def endpoint_rows(modulus, powers, m, target):
    """Assign prime powers to three unlabelled groups, directly lifting normalized parts."""
    def allocate(groups, index):
        if index == len(powers):
            if len(groups) != 3 or min(groups) <= m:
                return
            parts = []
            for group in groups:
                other = modulus // group
                part = other * pow(m*other, -1, group)
                if part >= target:
                    return
                parts.append(part)
            if sum(parts) == target:
                yield tuple(m*x for x in parts)
            return
        w = powers[index]
        for i in range(len(groups)):
            yield from allocate(groups[:i] + [groups[i]*w] + groups[i+1:], index + 1)
        if len(groups) < 3:
            yield from allocate(groups + [w], index + 1)
    yield from allocate([], 0)


def exclude_endpoint(n, power, counts):
    shift = n - 1
    if power == 1 or 27*shift**2 > power**3:
        return
    counts["endpoint weights to factor"] += 1
    m = n // power
    factors = prime_powers(shift)
    if len(factors) < 3 or max(factors.values()) > shift // (m+1)**2:
        return
    for row in endpoint_rows(shift, tuple(factors.values()), m, power):
        counts["endpoint CRT rows"] += 1
        if all(avoids(n, row, p) for p in factors):
            raise ArithmeticError(("shifted primes did not eliminate row", n, row))


def required_product(n, a, b):
    result = 1
    for h in range(6):
        value = n - h
        removed = gcd(value, 60)
        for i in range(h + 1):
            partial = gcd(value, a - i)
            if partial > 1:
                for j in range(h + 1 - i):
                    removed = lcm(removed, gcd(partial, b-j))
        result *= (value // removed)**(6-h)
    return result


def exclude_without_endpoint(n, row, counts):
    a, b, c = row
    if a + b <= 10:
        counts["binary row has small sum"] += 1
        return
    # One unknown row has no 1. The other can cost 64 times as much.
    forced = required_product(n, a, b)
    ordinary = n**9 // 3**9
    with_unit = (64*n**9) // 3**9
    if forced > ordinary * with_unit:
        counts["two ordinary row bound"] += 1
        return
    r = n if n % 2 else n - 1
    while r > c and not fmpz(r).is_prime():
        r -= 2
    if r > c:
        gap = n - r
        if gap <= 10:
            counts["near prime has small gap"] += 1
            return
        if forced > with_unit * gap**15:
            counts["near prime product bound"] += 1
            return
    raise ArithmeticError(("non-endpoint selection unresolved", n, row))


def main():
    cli = ArgumentParser(description=__doc__)
    cli.add_argument("bits", type=int, default=128, nargs="?")
    bits = cli.parse_args().bits
    if bits < 3:
        cli.error("use at least three bit positions")
    counts = Counter()
    for r, s, t in combinations(range(bits), 3):
        row = (1 << r, 1 << s, 1 << t)
        n = sum(row)
        exclude_endpoint(n, row[0], counts)
        exclude_without_endpoint(n, row, counts)
        counts["all weights"] += 1
    print(dict(counts))


if __name__ == "__main__":
    main()
