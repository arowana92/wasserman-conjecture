"""Certify 1+p^a(1+p^i+p^(2i))(1+p^k) for all primes and positive exponents."""

from collections import Counter
from functools import lru_cache
from itertools import combinations
from math import gcd, lcm, prod

import sympy as sp
from flint import fmpz, fmpz_poly

from check_polynomial_powers import coeffs, divisor_bound
from check_positive_root_weights import X, positive_after


T, Y = sp.symbols('T Y')
GRID = tuple((i, j) for i in range(3) for j in range(2))
MULTIPLICITIES = (3, 2, 2)


def centered(block, v, w):
    powers = [(i*v[0]+j*v[1], i*w[0]+j*w[1]) for i, j in block]
    lo, low = map(min, zip(*powers))
    return sum(T**(i-lo)*Y**(j-low) for i, j in powers)


def allocations():
    weight = sp.Poly(centered(GRID, (1, 0), (0, 1)), T, Y, domain=sp.ZZ)
    for z in GRID:
        rest = tuple(point for point in GRID if point != z)
        for x in combinations(rest, 2):
            y = tuple(point for point in rest if point not in x)
            for block in (x, y):
                factor = sp.Poly(centered(block, (1, 0), (0, 1)), T, Y, domain=sp.ZZ)
                quotient, remainder = weight.div(factor, auto=False)
                if remainder.is_zero:
                    assert quotient*factor == weight
                    break
            else:
                yield x, y, z


@lru_cache(maxsize=None)
def projection(first, second, v, w):
    f, g = centered(first, v, w), centered(second, v, w)
    if T not in f.free_symbols:
        bound = f
    elif T not in g.free_symbols:
        bound = g
    else:
        bound = sp.resultant(f, g, T)
    poly = sp.Poly(bound, Y, domain=sp.ZZ)
    assert not poly.is_zero
    terms = poly.as_dict()
    lo, hi = min(terms)[0], max(terms)[0]
    return tuple(int(terms.get((i,), 0)) for i in range(lo, hi+1))


@lru_cache(maxsize=None)
def roots(f):
    return tuple(int(z) for z, _ in fmpz_poly(list(f)).roots() if z > 0)


def height(v):
    return 4*v[0]+2*v[1]


def sectors(first, second, pairs):
    assert first[0]*second[1]-first[1]*second[0] == 1
    mid = tuple(i+j for i, j in zip(first, second))
    orientations = []
    for v, w in ((first, second), (second, first)):
        fs = [projection(f, g, v, w) for f, g in pairs]
        degree = sum(k*(len(h)-1) for k, h in zip(MULTIPLICITIES, fs))
        orientations.append((v, w, fs, degree))
    if any(height(mid) <= degree+1 for _, _, _, degree in orientations):
        yield from sectors(first, mid, pairs)
        yield from sectors(mid, second, pairs)
    else:
        yield from orientations


def binary_offset(x, y, v, w, fs):
    def specialized(block):
        f = sum(2**(i*w[0]+j*w[1])*X**(i*v[0]+j*v[1]) for i, j in block)
        while f[0] == 0:
            f = f//X
        return f
    u = fmpz_poly(list(divisor_bound(coeffs(specialized(x)), coeffs(specialized(y)))))
    while all(c % 2 == 0 for c in u.coeffs()):
        u = u//2
    a, c = [abs(sum(value*2**j for j, value in enumerate(h))) for h in fs[1:]]
    assert a and c
    i, z = 2**w[0]*X**v[0], 2**w[1]*X**v[1]
    m = (1+i+i**2)*(1+z)
    assert positive_after(729*m**4-u**6*a**4*c**4, 2)


def direct(x, y, i, k, base):
    power, block = base**i, base**k
    m = (1+power+power**2)*(1+block)
    xval, yval = [sum(power**l*block**j for l, j in part) for part in (x, y)]
    u = gcd(xval, yval)
    shared = gcd(u, base)
    while shared > 1:
        u //= shared
        shared = gcd(u, shared)
    s = gcd(m, lcm(xval, yval))
    assert 27*m*m > u**3*s**2


def check_ray(pairs, v, counts):
    def specialize(block):
        f = sum(X**(i*v[0]+j*v[1]) for i, j in block)
        while not f[0]:
            f = f//X
        return coeffs(f)

    u, a, c = [fmpz_poly(list(divisor_bound(specialize(f), specialize(g)))) for f, g in pairs]
    power, block = X**v[0], X**v[1]
    m = (1+power+power**2)*(1+block)
    difference = 729*m**4-u**6*a**4*c**4
    assert difference and difference[difference.degree()] > 0
    cutoff = 2
    while not positive_after(difference, cutoff):
        cutoff *= 2
    for base in range(2, cutoff):
        direct(*pairs[0], *v, base)
        counts['small integer ray evaluations'] += 1
    counts['largest translation'] = max(counts['largest translation'], cutoff)


def check():
    counts = Counter()
    root_rays = set()
    positive_roots = set()
    for x, y, z in allocations():
        counts['retained shapes'] += 1
        pairs = [(x, y), (x, tuple(sorted(y+(z,)))), (y, tuple(sorted(x+(z,))))]
        ray_checks, ceilings = set(), {}
        for v, w, fs, degree in sectors((1, 0), (0, 1), pairs):
            counts['oriented sectors'] += 1
            for value in set(r for f in fs for r in roots(f)):
                positive_roots.add(value)
                if value == 1:
                    if all(v):
                        ray_checks.add(v)
                        root_rays.add(v)
                else:
                    assert value == 2
                    binary_offset(x, y, v, w, fs)
                    counts['binary offset checks'] += 1
            norm = prod(sum(abs(c) for c in f)**k for f, k in zip(fs, MULTIPLICITIES))
            gap = height(v)+height(w)-degree
            assert gap > 0
            b = 0
            while 27*2**(b*height(v)) <= norm:
                ell = 0
                while 27*2**(b*height(v)+ell*gap) <= norm:
                    i = (b+ell)*v[0]+ell*w[0]
                    k = (b+ell)*v[1]+ell*w[1]
                    if i and k:
                        limit = norm if ell else prod(abs(sum(f))**j for f, j in zip(fs, MULTIPLICITIES))
                        exponent = b*height(v)+ell*gap
                        high = int(fmpz(limit//27).root(exponent))
                        ceilings[i, k] = max(ceilings.get((i, k), 0), high)
                    ell += 1
                b += 1
        for (i, k), high in ceilings.items():
            if high > 64:
                common = gcd(i, k)
                ray_checks.add((i//common, k//common))
        for v in sorted(ray_checks):
            check_ray(pairs, v, counts)
            counts['univariate shape and ray pairs'] += 1
        for (i, k), high in ceilings.items():
            common = gcd(i, k)
            if (i//common, k//common) not in ray_checks:
                assert high <= 64
                for p in sp.primerange(2, high+1):
                    direct(x, y, i, k, int(p))
                    counts['finite prime and row evaluations'] += 1
    assert counts['retained shapes'] == 42
    print(dict(counts))
    print('Positive integer roots', sorted(positive_roots))
    print('Interior rays with zero projections at 1', sorted(root_rays))


if __name__ == '__main__':
    check()
