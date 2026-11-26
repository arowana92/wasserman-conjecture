"""Exact checks of sparse binomial gcds and their specialization on exponent lines."""
from itertools import combinations_with_replacement, product
from math import gcd
from random import Random

from flint import fmpz_poly

X = fmpz_poly([0, 1])


def primitive(poly):
    content = gcd(*(abs(int(c)) for c in poly.coeffs()), 0)
    return poly // content


def one_line(p, d, positions, sign, slope, slopes):
    original = (X**d + sign).gcd(sum(X**a for a in positions))
    integer_gcd = gcd(p**d + sign, sum(p**a for a in positions))
    divisor = abs(int(original(p)))
    assert integer_gcd % divisor == 0

    first = p**d * X**slope + sign
    second = sum(p**a * X**v for a, v in zip(positions, slopes))
    common = primitive(first.gcd(second))
    specialized = abs(int(common(1)))
    assert divisor % specialized == 0
    return specialized, divisor


def line_checks():
    rng = Random(785)
    counts = [0, 0, 0]
    for p, s in ((5, 3), (5, 4), (7, 4), (7, 6), (11, 8), (17, 10)):
        for sign in (-1, 1):
            for _ in range(400):
                d = rng.randint(1, 12)
                positions = (0, *sorted(rng.randrange(13) for _ in range(s-1)))
                for slope, slopes in ((d, positions), (rng.randrange(7),
                                      (0, *(rng.randrange(7) for _ in range(s-1))))):
                    line, full = one_line(p, d, positions, sign, slope, slopes)
                    counts[0] += 1
                    counts[1] += line > 1
                    counts[2] += line < full

    # Opposite pairs have proportional gaps within each pair. Translate
    # one pair in slope while keeping both pairs at their original heights.
    for p, w, gap in product((5, 7), range(1, 5), range(1, 7)):
        line, full = one_line(p, 3*w, (0, w, gap, gap+w),
                              1, 3, (0, 1, 8, 9))
        assert line == p**w+1 and full % line == 0
    # Both matchings occur originally. Only the gap-5 matching survives
    # on the weighted line. Its value is a proper nonconstant divisor.
    for p in (5, 7, 11):
        line, full = one_line(p, 15, (0, 3, 5, 8), 1, 3, (0, 4, 1, 5))
        assert 1 < line == p**5+1 < full

    assert min(counts) > 0
    print(f"{counts[0]} lines, {counts[1]} nonconstant evaluated factors, "
          f"{counts[2]} proper divisors, 51 additional paired lines")


def quadratic_characters():
    for p in (3, 5, 7, 11, 17, 23):
        chi = [0] + [1 if pow(j, (p-1)//2, p) == 1 else -1 for j in range(1, p)]
        tau = fmpz_poly(chi)
        conjugate = fmpz_poly([0, *reversed(chi[1:])])
        phi = fmpz_poly([1]*p)
        assert (tau*conjugate) % phi == p
        assert all(sum(abs(a+b) for a in chi) >= p-1 for b in range(-p, p+1))
    checks = 0
    for p in (3, 5, 7, 11, 13):
        for t in range((p-1)//2):
            seen = set()
            for exponents in combinations_with_replacement(range(p), t):
                moments = tuple(sum(pow(a, k, p) for a in exponents) % p
                                for k in range(1, t+1))
                assert moments not in seen
                seen.add(moments)
                checks += 1
    print(f"Gauss sums at six primes and {checks} signed-moment multiset checks")


if __name__ == '__main__':
    line_checks()
    quadratic_characters()
