"""Verify the two-exponent family at all odd prime bases by exact elimination."""

from collections import Counter
from functools import lru_cache
from itertools import combinations
from math import gcd

import sympy as sp
from flint import fmpz, fmpz_poly

from check_polynomial_powers import coeffs, common_multiple, divisor_bound, trim
from check_positive_root_weights import X, numerical_test, positive_after, two_block_bounds, two_blocks


x, y = sp.symbols("x y")
B = x**2 + 2*x + 2
C = 2*x*(x+1)


def weight(ray):
    a, j = ray
    return 4*a + 2*j - 3*min(a, j)


def add(v, w):
    return v[0]+w[0], v[1]+w[1]


def determinant(v, w):
    return v[0]*w[1] - v[1]*w[0]


def subdivide(v, w, degree):
    """Return left endpoints of unimodular cones with enough growth."""
    assert determinant(v, w) == 1
    middle = add(v, w)
    if weight(middle) >= degree+2:
        return [v]
    return subdivide(v, middle, degree) + subdivide(middle, w, degree)


@lru_cache(maxsize=None)
def projection(poly, delta, v):
    n, m = v
    if n == 0:
        return poly
    content = abs(gcd(*poly))
    if len(poly) == 1:
        return (content,)
    primitive = sum((a//content)*x**i for i, a in enumerate(poly))
    relation = y*x**m*B**n - (delta-C)**n
    bound = sp.Poly(content*sp.resultant(primitive, relation, x), y)
    assert not bound.is_zero and bound.degree() < len(poly)
    # Retain the full homogeneous degree for the reverse projection.
    return tuple(int(bound.nth(i)) for i in range(len(poly)))


def directed(poly, delta, v, w):
    result = projection(poly, delta, v)
    if v[0] and determinant(v, w) < 0:
        result = result[::-1]
    return trim(result)


@lru_cache(maxsize=None)
def positive_integer_roots(poly):
    return {int(z) for z, _ in fmpz_poly(list(poly)).roots() if z > 0}


def finite_domain(v, w, u, s):
    degree = 3*(len(u)-1) + 2*(len(s)-1)
    gap = weight(add(v, w)) - degree
    assert gap > 0
    norm = sum(map(abs, u))**3 * sum(map(abs, s))**2
    ell = 0
    while 27*3**(gap*ell) <= norm:
        # On the ray (ell=0), use the exact values at 1.
        ceiling = norm if ell else abs(sum(u))**3 * abs(sum(s))**2
        b = 0
        while 27*3**(gap*ell + weight(v)*b) <= ceiling:
            a = (b+ell)*v[0] + ell*w[0]
            j = (b+ell)*v[1] + ell*w[1]
            if a and j:
                exponent = gap*ell + weight(v)*b
                yield (a, j), int(fmpz(ceiling//27).root(exponent))
            b += 1
        ell += 1


def strip_x(poly):
    assert poly
    while poly[0] == 0:
        poly = poly//X
    return poly


def check_ray(lower, upper, v, counts):
    """Univariate size inequality for every odd prime power on this ray."""
    a, j = v
    assert a > 0 and j > 0 and gcd(a, j) == 1
    parts = [h(X**a)*X**j + l(X**a) for l, h in zip(lower, upper)]
    n = sum(parts)
    m = n-1
    parts_free = [coeffs(strip_x(f)) for f in parts]
    u_free = divisor_bound(coeffs(strip_x(n)), *parts_free)
    u = X**min(v)*fmpz_poly(list(u_free))
    bins = tuple(divisor_bound(coeffs(m), i, k) for i, k in combinations(parts_free, 2))
    s = fmpz_poly(list(common_multiple(bins)))
    difference = 729*m**4 - u**6*s**4
    assert difference and difference[difference.degree()] > 0, (v, parts)
    cutoff = 3
    while not positive_after(difference, cutoff):
        cutoff *= 2
    for t in range(3, cutoff):
        numerical_test(tuple(int(f(t)) for f in parts))
        counts['small integer ray evaluations'] += 1
    counts['largest ray cutoff'] = max(counts['largest ray cutoff'], cutoff)


def check():
    counts = Counter()
    all_roots, root_rays = set(), set()
    for lower, upper in two_blocks():
        u_total, s_total = two_block_bounds(lower, upper)
        u, s = coeffs(u_total//(2*X)), coeffs(s_total)
        degree = 3*(len(u)-1) + 2*(len(s)-1)
        rays = subdivide((1, 0), (1, 1), degree) + subdivide((1, 1), (0, 1), degree) + [(0, 1)]
        ceilings = {}
        special = set()
        for left, right in zip(rays, rays[1:]):
            for v, w in ((left, right), (right, left)):
                first, second = directed(u, 0, v, w), directed(s, 1, v, w)
                roots = positive_integer_roots(first) | positive_integer_roots(second)
                all_roots.update(roots)
                # A positive power of an odd prime cannot equal 4 or 16.
                assert roots <= {1, 4, 16}
                if 1 in roots and all(v):
                    special.add(v)
                    root_rays.add(v)
                    counts['interior root incidences'] += 1
                for pair, high in finite_domain(v, w, first, second):
                    ceilings[pair] = max(ceilings.get(pair, 0), high)
                counts['oriented cones'] += 1

        # Large finite prime ranges are cheaper to check on an entire ray.
        for (a, j), high in ceilings.items():
            if high > 1000:
                g = gcd(a, j)
                special.add((a//g, j//g))
        for v in sorted(special):
            check_ray(lower, upper, v, counts)
            counts['univariate shape and ray pairs'] += 1
        for (a, j), high in ceilings.items():
            g = gcd(a, j)
            if (a//g, j//g) in special:
                continue
            assert high <= 1000
            for p in range(3, high+1, 2):
                if fmpz(p).is_prime():
                    power, block = p**a, p**j
                    numerical_test(tuple(int(h(power))*block + int(l(power)) for l, h in zip(lower, upper)))
                    counts['direct prime and row evaluations'] += 1
        counts['shapes'] += 1

    assert counts['shapes'] == 575
    print(dict(counts))
    print('Positive integer roots', sorted(all_roots))
    print('Interior rays with a zero at 1', sorted(root_rays))


if __name__ == '__main__':
    check()
