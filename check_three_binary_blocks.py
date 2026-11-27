"""Exact checks for the estimates for three binary blocks with two free gaps."""

from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations, product
from math import factorial, gcd, prod

from check_binary_gap_bounds import entropy_below


def minors(first, second, last):
    for i, j in combinations(range(3), 2):
        det = first[i] * second[j] - first[j] * second[i]
        if det:
            left = second[i] * last[j] - second[j] * last[i]
            right = first[j] * last[i] - first[i] * last[j]
            return det, left, right
    return None


def slope(a, b):
    for u in range(1, 33):
        v = (2 * u * b + a) // (2 * a)
        r = v * a - u * b
        if 32 * abs(r) <= a:
            return u, v, r
    raise ArithmeticError((a, b))


def odd_rational(t):
    t = abs(t)
    v = 0
    while t.numerator % 2 == 0:
        t /= 2
        v += 1
    while t.denominator % 2 == 0:
        t *= 2
        v -= 1
    return v, t


def check_matrix():
    count = 0
    # The signs in Cramer's rule and the three distinct targets e_i matter.
    for first in ((1, 0, 0), (3, 4, 0), (1, 4, 10)):
        for second in ((3, 4, 8), (0, 1, 0), (5, 0, 2)):
            for lower in ((1, 2, 4), (4, 1, 0), (0, 0, 7)):
                for slot in (-1, 0, 1, 2):
                    last = tuple(z - int(i == slot) for i, z in enumerate(lower))
                    data = minors(first, second, last)
                    if data is None:
                        continue
                    det, left, right = data
                    for a, b in ((5, 2), (13, 5), (8, 8), (15, 0)):
                        vals = [x * 2**a + y * 2**b + z for x, y, z in zip(first, second, last)]
                        w = gcd(*vals)
                        assert (det * 2**a - left) % w == (det * 2**b - right) % w == 0
                        count += 1
    for a in range(1, 381):
        for b in range(a + 1):
            u, v, r = slope(a, b)
            assert 0 <= v <= u <= 32 and abs(r) * 32 <= a
    first, second, last = (1, 128, 32768), (1, 256, 16384), (18, 1152, 737280)
    det, left, right = minors(first, second, last)
    x, y = Fraction(left, det), Fraction(right, det)
    assert (x, y) == (-27, 9) and x**2 == y**3
    assert all(x * f + y * g + t == 0 for f, g, t in zip(first, second, last))
    for k in range(24, 129):
        parts = [f * 2**(3*k) + g * 2**(2*k) + t for f, g, t in zip(first, second, last)]
        assert gcd(*parts) % (2**k + 3) == 0
        assert not any(z & w for z, w in combinations(parts, 2))
    return count


def check_roots():
    successes = 0
    for u in range(2, 7):
        for v in range(1, u):
            if gcd(u, v) > 1:
                continue
            for t in (Fraction(1), Fraction(1, 3), Fraction(3, 5), Fraction(5, 3)):
                for j, k in product(range(-2, 3), repeat=2):
                    x = -(Fraction(2)**j) * t**u
                    y = Fraction(2)**k * t**v
                    # Build positive integer coefficients from a rational root.
                    d = 2**4 * t.denominator**u
                    m = d
                    c = -d*x - m*y
                    if c.denominator != 1 or c <= 0:
                        continue
                    assert d % t.denominator == int(c) % t.numerator == 0
                    assert odd_rational(x) == (j, t**u)
                    assert odd_rational(y) == (k, t**v)
                    successes += 1
    return successes


def matrix_rank(differences, dimension):
    pivots = {}
    for vector in differences:
        v = list(map(Fraction, vector))
        for i in range(dimension):
            if not v[i]:
                continue
            if i not in pivots:
                pivots[i] = [z / v[i] for z in v]
                break
            multiple = v[i]
            v = [z - multiple * q for z, q in zip(v, pivots[i])]
        if len(pivots) == dimension:
            break
    return len(pivots)


def check_lattices():
    checks = Counter()
    for widths in ((4,), (4, 8), (2, 8), (2, 4, 8), (1, 4, 8), (4, 4, 4)):
        dim = len(widths)
        coeff = [2**(2*j+1) for j in range(dim)]
        for modulus in range(3, 46, 2):
            def label(vector):
                return sum(z * a for z, a in zip(vector, coeff)) % modulus
            points = defaultdict(list)
            for vector in product(*(range(w) for w in widths)):
                points[label(vector)].append(vector)
            # Points with distance greater than one cannot affect separation in this box.
            scale = max(widths)
            length = min((max(abs(z) * scale // w for z, w in zip(v, widths))
                          for v in product(*(range(-w, w+1) for w in widths))
                          if any(v) and label(v) == 0), default=2*scale)
            for group in points.values():
                rank = matrix_rank(([q-t for q, t in zip(p, group[0])] for p in group), dim)
                if rank == dim:
                    assert len(group) * modulus <= 3**dim * factorial(dim) * prod(widths)
                    checks["full rank"] += 1
                else:
                    assert len(group) * length**rank <= (length + 4 * scale)**rank
                    checks["smaller rank"] += 1
    return checks


if __name__ == "__main__":
    assert entropy_below(3, 5, 39, 40)
    assert entropy_below(1, 5, 29, 40)
    print("Matrix congruences:", check_matrix())
    print("Rational-root examples:", check_roots())
    print("Lattice translates:", check_lattices())
