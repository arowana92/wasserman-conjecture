"""Exact size comparisons and polynomial controls for the improved seven-unit cutoffs."""

from collections import Counter
from fractions import Fraction
from itertools import combinations, combinations_with_replacement
from math import gcd, lcm, prod

from check_dense_offsets import free_rows
from check_gap_resultants import resultant
from check_seven_large_bases import X, actual_row_checks, poly


def below_one(terms):
    """Compare a product of positive rational powers with one by clearing exponents."""
    power = lcm(*(denominator for _, _, denominator in terms))
    return prod(Fraction(base)**(numerator*(power//denominator))
                for base, numerator, denominator in terms) < 1


def numerical_bounds():
    b = 64
    for terms in (
        [(3, 3, 2), (9, 1, 1), (6, 2, 3), (b+1, 2, 3)],
        [(Fraction(19, 4), 4, 5), (9, 1, 1), (4, 8, 15), (b+1, 14, 15)],
        [(5, 3, 2), (16, 1, 1), (10, 8, 9), (b+1, 2, 9)],
        [(129, 3, 21), (16, 1, 1), (10, 2, 3), (b+1, 2, 3)],
        [(3, 3, 2), (54, 1, 1), (Fraction(b**6, 27), 1, 9)],
    ):
        assert below_one(terms + [(b, -2, 1)])
    small = 0
    for height in range(1, 21):
        for d in range(1, height+1):
            for k in range(1, d//3+1):
                if d%k or (d//k)%2 == 0:
                    continue
                u6 = min(5**(3*d), (2**(d//k)+1)**6)
                ac2 = 10**(d-k)*(b+1)**(2*k)*16**height
                assert u6*ac2**2 < (27*b**(2*height))**2
                small += 1
    b = 23
    for terms in (
        [(3, 1, 1), (6, 1, 3), (b+1, 1, 3)],
        [(3, 1, 1), (2, 8, 15), (b+1, 7, 15)],
        [(4, 1, 1), (10, 1, 2)],
        [(4, 1, 1), (6, 1, 3), (b+1, 1, 3)],
    ):
        assert below_one(terms + [(2, 1, 21), (b, -1, 1)])
    print(f'Cutoff inequalities and {small} small-height triples.')


def four_positions():
    counts = Counter()
    for d in range(1, 26):
        f = X**d+1
        for offsets in combinations_with_replacement(range(24), 3):
            places = (0, *offsets)
            if len(set(places)) == 1:
                continue
            h = f.gcd(poly(places))
            k = h.degree()
            if k == d:
                continue
            distinct = len({n%d for n in places}) == 4
            if not distinct and k:
                assert 3*k <= d and h == X**k+1
                counts['proper collision factors'] += 1
            for z in set(places):
                y = list(places)
                y.remove(z)
                residues = {n%d for n in y}
                if len(residues) == 2 and k:
                    pair = next((i, j) for i, j in combinations(range(3), 2)
                                if y[i]%d == y[j]%d)
                    i, j = pair
                    if (y[i]-y[j])%(2*d):
                        continue
                    other = y[3-i-j]
                    assert (other-y[i])%k == 0
                    reduced = (poly(places) % f).coeffs()
                    content = gcd(*(int(c) for c in reduced))
                    assert sum((int(c)//content)**2 for c in reduced) <= 6
                    determinant = resultant(f, poly(tuple(y)))
                    assert determinant > 0
                    for b in (3, 5, 17, 64):
                        u = gcd(b**d+1, sum(b**e for e in y))
                        assert u <= 2**(d//k)+1
                    counts['two-residue scale bounds'] += 1
                if len(residues) == 3 and 3*k > d:
                    assert distinct and f.gcd(poly(tuple(y))) == 1
                    # At every root of the four-term gcd the value of y is -X^z.
                    assert not (poly(tuple(y))+X**z) % h
                    determinant = resultant(f, poly(tuple(y)))
                    assert determinant**2*(d-k)**(d-k) <= (3*d-k)**(d-k)
                    counts['large common factors'] += 1
    assert all(counts[name] > 0 for name in
               ('proper collision factors', 'two-residue scale bounds', 'large common factors'))
    print('Four-position checks:', dict(counts))


def binomial_resultants():
    counts = Counter()
    for d in range(1, 33):
        for ell in range(1, 33):
            height = max(d, ell)
            for first in (X**d+2, 2*X**d+1):
                for second in (X**ell+2, 2*X**ell+1):
                    if first == second:
                        continue
                    r = resultant(first, second)
                    g = gcd(d, ell)
                    # Even resultant content and multiplicity are irrelevant at odd C.
                    b1, b2 = (1 if f[0] == 2 else -1 for f in (first, second))
                    u, v = b1*(ell//g), b2*(d//g)
                    delta = abs((-1)**(ell//g)*2**(u-min(u,v)) -
                                (-1)**(d//g)*2**(v-min(u,v)))
                    assert 0 < delta <= 4**height
                    assert (r//(r & -r)) == (delta//(delta & -delta))**g
                    counts['binomial pairs'] += 1
            for s in range(1, ell):
                if len({0, s%d, ell%d}) != 3:
                    continue
                for f in (X**d+2, 2*X**d+1):
                    r = resultant(f, 1+X**s+X**ell)
                    assert 0 < r <= 3**height
                    counts['weighted trinomial resultants'] += 1
    print('Resultant checks:', dict(counts))


def mixed_rows():
    counts = Counter()
    for p, step in ((3, 3), (5, 2), (23, 1), (29, 1)):
        for height in (21, 30, 42, 63):
            for middle in combinations_with_replacement((0, 1, 3, 5, height//3, height-1), 4):
                levels = (0, *middle, height)
                if any(levels.count(v) >= p for v in set(levels)):
                    continue
                m = sum(p**(step*h) for h in levels)
                for shift in (1, 2):
                    n = 1+p**shift*m
                    if n%4 != 3:
                        continue
                    counts['mixed weights'] += 1
                    for row in free_rows(n, p):
                        counts['mixed rows'] += 1
                        if gcd(*row) < 7:
                            continue
                        counts['mixed scales at least seven'] += 1
                        bins = prod(gcd(m, v, w) for v, w in combinations(row, 2))
                        assert (m//2) % (bins//(bins & -bins)) == 0
                        assert bins//(bins & -bins) < m//2
    assert counts['mixed scales at least seven'] > 0
    print('Mixed cutoff comparisons:', dict(counts))


if __name__ == '__main__':
    numerical_bounds()
    binomial_resultants()
    four_positions()
    mixed_rows()
    actual_row_checks(((3, 4), (5, 3), (67, 1), (2, 6)))
