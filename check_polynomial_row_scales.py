"""Check the polynomial recovery rule against integer gcds and direct divisibility."""

from functools import lru_cache
from itertools import combinations
from math import gcd, lcm
from random import Random

from sympy import Poly, QQ, ZZ

from check_polynomial_row_divisors import X, bezout_bound
from check_primitive_row_scales import profile


def poly(expression):
    return Poly(expression, X, domain=QQ)


@lru_cache(maxsize=None)
def polynomial_gcd_data(polys):
    denominator = lcm(*(int(v.q) for f in polys for v in f.all_coeffs()))
    integral = tuple((f * denominator).set_domain(ZZ) for f in polys)
    common, bound = bezout_bound(integral)
    return common.set_domain(QQ) / denominator, bound, tuple(f.exquo(common) for f in integral)


def profile_on_class(vectors, t):
    period = 1

    def gcd_piece(polys):
        nonlocal period
        common, bound, quotients = polynomial_gcd_data(tuple(polys))
        period = lcm(period, bound)
        multiplier = gcd(bound, *(int(f.eval(t)) for f in quotients))
        return poly(common * multiplier)

    rows = []
    for row in vectors:
        parts = tuple(map(poly, row))
        common = gcd_piece(parts)
        rows.append(tuple(f.exquo(common) for f in parts))
    u, v = (sum(row) for row in rows)
    base_gcd = gcd_piece((u, v))
    h, j = (gcd_piece(tuple(a*b for a, b in combinations(row, 2))) for row in rows)
    H = (h*j).exquo(gcd_piece((h, j)))
    return period, u*v, H, base_gcd


def branch(n, H):
    quotient, remainder = H.div(n)
    if remainder.is_zero or remainder.LC() >= 0:
        return None
    scale, tail = quotient.div(-remainder)
    if not tail.is_zero or scale.is_zero or scale.LC() <= 0:
        return None
    return scale, -remainder


def predicted(n, H, t):
    result = branch(n, H)
    if result:
        scale, complement = (f.eval(t) for f in result)
        if scale.q == complement.q == 1 and min(scale, complement) > 0:
            return int(scale)
    return None


def scalar_controls():
    rng = Random(283)
    n = poly(127*X**4+97*X+53)
    cases = [(poly(a)*X**2+poly(b)*X+poly(c), poly(r))
             for a in (0, 1) for b in (1, 3) for c in (0, 2)
             for r in (0, 1, -1, -2, -X, -X-1, -X**3)]
    for _ in range(36):
        d = rng.randrange(1, 8)
        k = poly((X+rng.randrange(5))/d)
        g = poly((X+rng.randrange(5))/rng.randrange(1, 6))
        cases.extend([(k*g, -k), (k*g+poly(1), -k)])
    tests = candidates = 0
    for q, r in cases:
        H = n*q+r
        for t in range(100000, 100060):
            value = H.eval(t)
            if value.q != 1:
                continue
            m, value = int(n.eval(t)), int(value)
            k = -value % m
            g = None
            if k and value % k == 0:
                weight = 1+value//k
                if weight % m == 0:
                    g = weight//m
            assert g == predicted(n, H, t)
            tests += 1
            candidates += g is not None
    return tests, candidates


def row_controls():
    h = 7*X
    d = h+1
    e = h*d+1
    f = h*d*e-1
    vectors = [
        ((1, 2, 4), (e*f, d*f, d*e)),
        ((X+2, 2*X+5, 5*X+1), (3*X+1, 9*X+3, 7*X+5)),
        ((2*X+2, 6*X+6, 3*X+1), (3*X-1, 6*X-2, 5*X+4)),
        ((X*(X+1), X*(X+2), X*(X+5)), ((X+3)**2, X*(X+3), X**2)),
        ((X+11, X+23, X+35), (4*X+1, 8*X+11, 12*X+21)),
    ]
    evaluations = carries = 0
    for i, pair in enumerate(vectors):
        rows = [tuple(map(poly, row)) for row in pair]
        for anchor in range(1, 19):
            period, n, H, c = profile_on_class(pair, anchor)
            for j in (100, 999, 10000):
                t = anchor+j*period
                directions = []
                for row in rows:
                    values = tuple(int(z.eval(t)) for z in row)
                    directions.append(tuple(a//gcd(*values) for a in values))
                u, v = map(sum, directions)
                a, b = (gcd(*(z*w for z, w in combinations(row, 2))) for row in directions)
                assert (n.eval(t), H.eval(t), c.eval(t)) == (u*v, lcm(a, b), gcd(u, v))
                small, _, weight = profile(*directions)
                actual = weight//small if weight else None
                assert actual == (predicted(n, H, t) if gcd(u, v) == 1 else None)
                if i == 0:
                    assert actual == t
                    for g in (1, 7, t, 7**4*t, 123456789):
                        for row in directions:
                            scale = g*small//sum(row)
                            parts = [scale*a for a in row]
                            N = sum(parts)
                            valuation = 0
                            while N:
                                N //= 7
                                parts = [a//7 for a in parts]
                                valuation += N-sum(parts)
                            assert valuation >= 1
                            carries += 1
                evaluations += 1
    return evaluations, carries


if __name__ == "__main__":
    print("Scalar evaluations and surviving scales:", scalar_controls())
    print("Primitive direction profiles and factorial-floor controls:", row_controls())
