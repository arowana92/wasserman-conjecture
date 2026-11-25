"""Test odd weights with six binary ones, separating the endpoint exactly."""

from argparse import ArgumentParser
from collections import Counter
from itertools import combinations
from math import comb, gcd, prod

from check_sparse_odd import gap_suffices, partitions
from check_three_bit_range import exclude_without_endpoint
from search_endpoint import find_row
from search_even_endpoint import avoids, prime_powers


def endpoint_shapes(units):
    m = sum(units)
    for i, c in enumerate(units):
        rest = units[:i] + units[i+1:]
        for j in (1, 2, 3):
            a = rest[0] + rest[j]
            b = m - c - a
            upper_u = gcd(a, b)
            upper_u //= upper_u & -upper_u
            upper_s = gcd(a*b, m)
            if upper_u >= 3 and 27*(m//upper_s)**2 <= upper_u**3:
                yield a, b, c, upper_u, upper_s


def endpoint(n, q, m, shapes, counts):
    for a, b, c, upper_u, upper_s in shapes:
        unit_part = q*c + 1
        u = gcd(upper_u, unit_part)
        s = gcd(upper_s, unit_part)
        if u < 3 or 27*(m//s)**2 > u**3:
            continue
        counts["binary rows requiring factorization"] += 1
        row = (q*a, q*b, unit_part)
        full = [(p, w) for p, w in prime_powers(n).items() if not avoids(n, row, p)]
        shifted = [(p, w) for p, w in prime_powers(m).items() if not avoids(n, row, p)]
        scale = prod(w for _, w in full)
        shift = prod(w for _, w in shifted)
        other = find_row(n, scale, full, shift, shifted)
        if other:
            raise ArithmeticError(("potential counterexample", n, row, other))


def check(bits):
    counts = Counter()
    for exponents in combinations(range(1, bits-1), 4):
        units = (1,) + tuple(1 << j for j in exponents)
        m = sum(units)
        shapes = list(endpoint_shapes(units))
        for k in range(1, bits-m.bit_length()+1):
            q = 1 << k
            n = m*q + 1
            counts["weights"] += 1
            if n <= 1_000_000:
                counts["previous million-weight theorem"] += 1
                continue
            endpoint(n, q, m, shapes, counts)
            if gap_suffices(n):
                counts["gap and height six"] += 1
                continue
            counts["weights requiring digit partitions"] += 1
            actual = (1,) + tuple(q*x for x in units)
            for row in partitions(actual):
                if len(row) == 3:
                    exclude_without_endpoint(n, sorted(row), counts)
    assert counts["weights"] == comb(bits-1, 5)
    print(dict(counts))


if __name__ == "__main__":
    cli = ArgumentParser(description=__doc__)
    cli.add_argument("bits", type=int, nargs="?", default=64)
    check(cli.parse_args().bits)
