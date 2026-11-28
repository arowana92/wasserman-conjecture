"""Finite arithmetic checks for the unbounded absolute carry-depth theorem."""

from itertools import product
from math import gcd, isqrt


def prime(n):
    return n >= 2 and all(n % i for i in range(2, isqrt(n) + 1))


def first_carry(row, p):
    position = 0
    while any(row):
        position += 1
        if sum(x % p for x in row) >= p:
            return position
        row = tuple(x // p for x in row)
    return 0


def scalar_test(z, m, k, s):
    modulus = s
    while modulus <= m * z:
        if 2 * ((k * z) % modulus) > (m * z) % modulus:
            return False
        modulus *= s
    return True


def check_residues():
    checked = 0
    for s in (5, 7, 11, 13):
        for z in range(1, 160):
            for m in (3, 9, 27):
                for k in range(1, m // 2 + 1):
                    row = (k*z, k*z, (m-2*k)*z)
                    assert scalar_test(z, m, k, s) == (first_carry(row, s) == 0)
                    checked += 1
    return checked


def check_grid():
    # Two one-digit blocks and a final digit testing {9*x} >= 2/3.
    denominator = 27
    for dimension in (1, 2):
        good = 0
        for point in product(range(denominator), repeat=dimension):
            successful = [k for k in (1, 3)
                          if all((k*x) % denominator < 9 for x in point)]
            if successful and all((9*x) % denominator >= 18 for x in point):
                good += 1
                for k in successful:
                    assert all(2*((k*x) % denominator) <= (9*x) % denominator
                               for x in point)
        aligned_words = 3**(2*dimension) - (3**dimension - 1)**2
        assert good == aligned_words  # Each last digit must equal two.


def check_sieve_period():
    a, c, s, period = 4, 27, 11, 121
    selected = {t for t in range(period) if 3*t % s < 4 and 3*t % period < 60}
    for divisor in (1, 5, 7, 35):
        length = a * divisor * period
        count = sum((c*s*t-1) % a == 0 and t % period in selected
                    and t*(c*s*t-1) % divisor == 0 for t in range(length))
        roots = (2 if divisor % 5 == 0 else 1) * (2 if divisor % 7 == 0 else 1)
        assert count == len(selected) * roots


def check_example(h, m, s, t, k):
    a, scale = 2**(h+1), 3**h
    c, v = scale*m, s*t
    n = c*v
    q, remainder = divmod(n-1, a)
    assert remainder == 0 and all(prime(p) for p in (q, s, t))
    assert len({2, 3, q, s, t}) == 5 and q > a and t > c*s
    assert gcd(q, n) == 1 and 0 < 2*k < m and k % 2 == 1
    cut = k*scale*v
    rows = ((1, 1, n-2), ((n-1)//2, (n-1)//2, 1), (cut, cut, n-2*cut))
    assert all(min(row) > 0 and sum(row) == n for row in rows)
    assert first_carry(rows[2], s) == first_carry(rows[2], t) == 0
    assert scalar_test(scale*t, m, k, s)
    assert not scalar_test(scale*t, m, 1, s)
    depths = []
    for p in (2, 3, q, s, t):
        at_prime = [first_carry(row, p) for row in rows]
        if min(at_prime):
            depths.append(max(at_prime))
        for row, first in zip(rows, at_prime):
            power = p
            floors = []
            while power <= n:
                if n // power > sum(x // power for x in row):
                    floors.append(power)
                power *= p
            assert (p**first if first else None) == next(iter(floors), None)
    assert min(depths) == h + 1
    print('Composite example:', n, 'base', s, 'multiplier', k, 'first common digit', h+1)


if __name__ == '__main__':
    print('Scalar residue tests:', check_residues())
    check_grid()
    check_sieve_period()
    check_example(1, 9, 53, 4099, 3)
    check_example(2, 27, 251, 62233, 9)
