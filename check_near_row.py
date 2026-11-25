"""Exact finite checks for the four-prime and small-part arguments."""

from fractions import Fraction
from itertools import combinations
from math import gcd, isqrt

from check_endpoint_pair_bound import rows


def is_prime(n):
    return n > 1 and all(n % d for d in range(2, isqrt(n) + 1))


def four_prime_example():
    n = 305
    data = ((151, 2), (101, 3), (61, 5), (43, 7))
    h = {p: n - p*u for p, u in data}
    assert all(is_prime(p) and 0 <= h[p] <= 4 < p for p, _ in data)
    for (_, u), (_, v) in combinations(data, 2):
        assert gcd(u, v) == 1 and 2*4*u*v < n
    assert max(p for p in range(2, n+1) if is_prime(p)) == 293

    seen = set()
    count = 0
    for row, coefficient in rows(n):
        count += 1
        avoided = [p for p in h if coefficient % p]
        for p, q in combinations(avoided, 2):
            assert sum(row[:2]) <= min(h[p], h[q])
            seen.add(row)
    assert seen == {(1, 1, 303), (1, 2, 302)}
    print(f"Four-prime example: {count} rows, {len(seen)} avoid two primes.")


def parts_one_and_three():
    selections = 0
    for n in range(3, 101):
        table = list(rows(n))
        first = [k for row, k in table if 1 in row]
        second = [k for row, k in table if 3 in row]
        common = {gcd(gcd(a, b), n*(n-1)*(n-2)) for a in first for b in second}
        assert all(gcd(pair_gcd, k) > 1 for pair_gcd in common for _, k in table)
        selections += len(first)*len(second)*len(table)
    print(f"Parts 1 and 3: {selections} ordered selections through weight 100.")


def degree_thresholds():
    half = Fraction(1, 2)

    def f(z, count):
        return z / (2*(count+z))

    # Each expression increases strictly in z.
    for z, left, middle, constant in (
        (Fraction(6, 7), 2, 2, Fraction(1, 4)),
        (Fraction(7, 8), 2, 1, Fraction(1, 6)),
        (Fraction(1), 1, 2, Fraction(1, 6)),
    ):
        assert f(half, left) + f(z, middle) + constant == half
    print("Exact cutoff thresholds: 6/7, 7/8, 1.")


if __name__ == "__main__":
    four_prime_example()
    parts_one_and_three()
    degree_thresholds()
