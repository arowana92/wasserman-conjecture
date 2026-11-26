"""Exact checks for the common-base resultant argument and its vanishing cases."""

from itertools import product
from math import gcd, prod

import sympy as sp


Y = sp.symbols("X Y Z")
T = sp.Symbol("T")


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def directional(poly, v):
    terms = sp.Poly(poly, *Y).terms()
    low = min(dot(index, v) for index, c in terms)
    return sum(c*prod(y**i for y, i in zip(Y, index))*T**(dot(index, v)-low)
               for index, c in terms)


def project(f, g, v):
    first, second = directional(f, v), directional(g, v)
    if T not in first.free_symbols:
        result = first
    elif T not in second.free_symbols:
        result = second
    else:
        result = sp.resultant(first, second, T)
    assert result != 0
    terms = sp.Poly(result, *Y).terms()
    assert len({dot(index, v) for index, c in terms}) == 1
    return [(index, int(c)) for index, c in terms]


def evaluate_mod_p_units(terms, p, exponents):
    powers = [(dot(index, exponents), int(c)) for index, c in terms]
    low = min(e for e, c in powers)
    value = sum(c*p**(e-low) for e, c in powers)
    return value, powers


def away_from(value, p):
    while value and value % p == 0:
        value //= p
    return abs(value)


def check():
    x, y, z = Y
    pairs = [(1+x, 1+y+z), (x+y-z, 2*x-z), (1+x*y, x*y+z-1)]
    vectors = [(3, 2, 4), (-2, 1, 0), (1, 1, 1), (0, 0, -1), (1, -1, 0)]
    direct, zero = 0, 0
    for f, g in pairs:
        f_terms, g_terms = [sp.Poly(poly, *Y).terms() for poly in (f, g)]
        for v in vectors:
            terms = project(f, g, v)
            norm = sum(abs(c) for index, c in terms)
            for e in product(range(-2, 4), repeat=3):
                for p in (2, 3, 5, 11):
                    f_value, _ = evaluate_mod_p_units(f_terms, p, e)
                    g_value, _ = evaluate_mod_p_units(g_terms, p, e)
                    bound, powers = evaluate_mod_p_units(terms, p, e)
                    common = away_from(gcd(f_value, g_value), p)
                    if common:
                        assert bound % common == 0
                    else:
                        assert bound == 0
                    ordered = sorted(exp for exp, c in powers)
                    if bound == 0:
                        gap = ordered[-1]-ordered[-2]
                        assert p**gap <= norm
                        zero += 1
                    else:
                        assert abs(bound) <= norm*p**(ordered[-1]-ordered[0])
                    direct += 1
    terms = project(1+x, 1+y+z, (3, 2, 4))
    r = sum(c*prod(variable**i for variable, i in zip(Y, index)) for index, c in terms)
    target = x**4+x**2*y**3+z**3-3*x**2*y*z
    assert sp.expand(r-target) == 0 or sp.expand(r+target) == 0
    for t in range(2, 1002):
        assert gcd(1+t**3, 1+t**2+t**4) == t**2-t+1
    print('Divisibility and coefficient bounds', direct)
    print('Vanishing projections', zero)
    print('Exact common factors on the ray (3,2,4)', 1000)


if __name__ == '__main__':
    check()
