"""Certify all primitive six-exponent patterns through normalized span 20.

Run with python-flint. All polynomial and integer calculations are exact.
"""

from collections import Counter
from functools import lru_cache
from itertools import combinations, combinations_with_replacement, product
from math import comb, gcd, lcm

from flint import fmpq_poly as RationalPoly, fmpz_poly as Poly

from check_dense_offsets import free_rows, has_size_obstruction

X = Poly([0, 1])
MAX_SPAN = 20


def odd_part(value):
    return value // (value & -value)


@lru_cache(maxsize=8192)
def polynomial(positions):
    return sum((X**i for i in positions), Poly(0))


def trim(positions):
    return polynomial(tuple(i - positions[0] for i in positions))


def exact_quotient(dividend, divisor):
    quotient, remainder = divmod(dividend, divisor)
    assert not remainder
    return quotient


@lru_cache(maxsize=8192)
def bezout_bound(first, second):
    f, k = Poly(list(first)), Poly(list(second))
    h = f.gcd(k)
    a = RationalPoly(exact_quotient(f, h))
    b = RationalPoly(exact_quotient(k, h))
    common, s, t = a.xgcd(b)
    c = lcm(int(s.denom()), int(t.denom()))
    assert common == 1 and s*a + t*b == 1
    assert (c*s).denom() == 1 and (c*t).denom() == 1
    return h, odd_part(c)


def bound(first, second):
    return bezout_bound(tuple(first.coeffs()), tuple(second.coeffs()))


def allocations(shape):
    seen = set()
    for pair in combinations(range(6), 2):
        x = tuple(shape[i] for i in pair)
        if x[0] == x[1]:
            continue
        remaining = tuple(shape[i] for i in range(6) if i not in pair)
        for size in (2, 3):
            for indices in combinations(range(4), size):
                y = tuple(remaining[i] for i in indices)
                if y[0] == y[-1] or (size == 2 and x > y):
                    continue
                z = tuple(remaining[i] for i in range(4) if i not in indices)
                allocation = x, y, z
                if allocation not in seen:
                    seen.add(allocation)
                    yield allocation


def shifted_positive(test):
    assert test.degree() >= 0 and test[test.degree()] > 0, test
    start = 3
    while test(start) <= 0 or min(test(X + start).coeffs()) < 0:
        start *= 2
    return start


def two_pairs(m, x, y):
    d, e = x[1]-x[0], y[1]-y[0]
    g = gcd(d, e)
    if (d//g) % 2 == 0 or (e//g) % 2 == 0 or m.degree() >= d+e-g:
        return 0
    common = X**g + 1
    q = exact_quotient(X**d + 1, common) * exact_quotient(X**e + 1, common)
    h, c = bound(m, q)
    start = shifted_positive(exact_quotient(m, h) - c*X**g)
    for b in range(3, start, 2):
        v = int(m(b))
        assert v >= b**g * gcd(v, int(q(b))), (m, x, y, b)
    return start


def pair_and_triple(m, x, y):
    f, k = trim(x), trim(y)
    u, cu = bound(f, k)
    shared = u.gcd(m)
    while shared.degree() > 0:
        u = exact_quotient(u, shared)
        shared = u.gcd(m)
    if u.degree() == 0 and cu*int(u[0]) < 7:
        return 0
    hx, cx = bound(m, f)
    hy, cy = bound(m, k)
    h = hx * exact_quotient(hy, hx.gcd(hy))
    reduced_m = exact_quotient(m, h)
    start = shifted_positive(27*reduced_m**2 - cu**3*(cx*cy)**2*u**3)
    for b in range(3, start, 2):
        v = int(m(b))
        fv, kv = int(f(b)), int(k(b))
        scale = odd_part(gcd(fv, kv))
        while gcd(scale, v) > 1:
            scale //= gcd(scale, v)
        s = odd_part(lcm(gcd(v, fv), gcd(v, kv)))
        assert scale < 7 or 27*v*v > scale**3*s*s, (m, x, y, b)
    return start


def count_shapes(span):
    # Independently count multisets from their occupied support and multiplicities.
    return sum(comb(5, k) for k in (3, 4, 5)
               for support in combinations(range(1, span+1), k) if gcd(*support) == 1)


def direct_controls():
    weights = set()
    for tail in combinations_with_replacement(range(5), 5):
        shape = (0,) + tail
        if len(set(shape)) < 4:
            continue
        expected = set()
        for labels in product(range(3), repeat=6):
            x, y, z = (tuple(e for e, lab in zip(shape, labels) if lab == group)
                       for group in range(3))
            if (len(x) == 2 and len(y) in (2, 3) and z and len(set(x)) >= 2
                    and len(set(y)) >= 2 and (len(y) == 3 or x <= y)):
                expected.add((x, y, z))
        assert set(allocations(shape)) == expected
        for p in (3, 5, 7):
            for shift, step in ((1, 1), (3, 2), (2, 5)):
                n = 1 + sum(p**(shift+step*e) for e in shape)
                weights.add((n, p))
    rows = 0
    for n, p in weights:
        for row in free_rows(n, p):
            assert gcd(*row) < 7 or has_size_obstruction(n, row), (n, p, row)
            rows += 1
    print(f'Direct controls: {len(weights)} weight/base pairs and {rows} carry-free rows.')


def scan(span):
    counts = Counter()
    largest_start = 0
    shapes = 0
    for tail in combinations_with_replacement(range(span+1), 5):
        shape = (0,) + tail
        if len(set(shape)) < 4 or gcd(*shape) != 1:
            continue
        shapes += 1
        m = polynomial(shape)
        for x, y, z in allocations(shape):
            counts['allocations'] += 1
            if len(y) == 2:
                start = two_pairs(m, x, y)
            else:
                start = pair_and_triple(m, x, y)
            if start:
                counts['two pairs' if len(y) == 2 else 'pair and triple'] += 1
                counts['integer tests'] += len(range(3, start, 2))
                largest_start = max(largest_start, start)
    assert shapes == count_shapes(span)
    print(f'Normalized span at most {span}: {shapes} shapes, {dict(counts)}, '
          f'largest polynomial threshold {largest_start}.')


if __name__ == '__main__':
    direct_controls()
    scan(MAX_SPAN)
