"""Verify polynomial and carry arguments for seven predecessor digit units."""

from collections import Counter
from functools import lru_cache
from itertools import combinations, combinations_with_replacement, product
from math import gcd

from flint import fmpz_poly as Poly

from check_dense_offsets import free_rows

X = Poly([0, 1])


@lru_cache(maxsize=8192)
def poly(positions):
    return sum((X**i for i in positions), Poly(0))


def bounds():
    b = 256
    assert 3**45 * 6**16 * 15**30 * (b+1)**28 < b**60
    assert 5**9 * 10**4 * 25**6 * (b+1)**4 < b**12
    assert 3**9 * (6*15)**6 < b**8
    assert 28*9**6 < 16**6 and 100*9**21 < 16**21 and 27*10**21 < 16**21


def common_trinomials():
    forms = [((0, a, b), poly((0, a, b))) for b in range(1, 37) for a in range(b+1)]
    powers = {d: X**(3*d)-1 for d in range(1, 37)}
    total = nonconstant = 0
    for i, (positions, p) in enumerate(forms):
        for _, q in forms[i+1:]:
            h = p.gcd(q)
            if h.degree() > 0:
                nonconstant += 1
                for a, b in combinations(positions, 2):
                    if b > a:
                        assert powers[b-a] % h == 0
            total += 1
    assert nonconstant
    print(f'{total} distinct trinomial pairs, {nonconstant} nonconstant gcds.')


def matching_factor(d, gaps):
    if not all(gaps):
        return Poly(1)
    j = gcd(d, *gaps)
    if any((e//j) % 2 == 0 for e in (d, *gaps)):
        return Poly(1)
    return X**j+1


def four_term_checks():
    counts = Counter()
    attained = False
    for d in range(1, 25):
        f = X**d+1
        for a, b, c in combinations_with_replacement(range(19), 3):
            terms = (0, a, b, c)
            actual = f.gcd(poly(terms))
            predicted = Poly(1)
            matches = 0
            for gaps in ((a, c-b), (b, c-a), (c, b-a)):
                h = matching_factor(d, gaps)
                matches += h.degree() > 0
                predicted *= h//h.gcd(predicted)
            assert matches <= 2 and actual == predicted
            counts['four-term gcds'] += 1
            if actual == f:
                continue
            k = actual.degree()
            assert 15*k <= 7*d
            attained |= 15*k == 7*d
            for shifted in set(terms):
                remaining = list(terms)
                remaining.remove(shifted)
                y = poly(tuple(remaining))
                u = f.gcd(y)
                if u.degree() > 0:
                    t = gcd(d, remaining[1]-remaining[0], remaining[2]-remaining[0])
                    assert u == X**(2*t)-X**t+1 and remaining[-1]-remaining[0] >= 4*t
                    assert actual == 1 and (d//t)%6 == 3
                    counts['cyclotomic row scales'] += 1
                residues = [i % (2*d) for i in remaining]
                if len(set(residues)) < 3:
                    assert 3*k <= d
                column = (y+X**shifted) % f
                assert sum(c*c for c in column.coeffs()) <= 16
                if len({i%d for i in remaining}) == 3:
                    assert sum(c*c for c in column.coeffs()) <= 6
                elif len({i%d for i in remaining}) == 2:
                    assert sum(c*c for c in column.coeffs()) <= 10
    assert attained and min(counts.values()) > 0
    print('Polynomial matching controls:', dict(counts))


def actual_row_checks():
    counts = Counter()
    for p, step in ((3, 6), (5, 4), (17, 2), (257, 1)):
        for tail in combinations_with_replacement(range(5), 5):
            shape = (0,) + tail
            if max(shape.count(e) for e in set(shape)) >= p:
                continue
            m = sum(p**(step*h) for h in shape)
            for a in (1, 2):
                n = 1 + p**a*m
                counts['weights'] += 1
                for row in free_rows(n, p):
                    counts['rows'] += 1
                    u = gcd(*row)
                    if u < 7:
                        continue
                    s = 1
                    for first, second in combinations(row, 2):
                        s *= gcd(first, second, m)
                    while s % 2 == 0:
                        s //= 2
                    assert 27*m*m > u**3*s*s, (n, p, row)
                    counts['large row gcds'] += 1
    assert min(counts.values()) > 0
    print('Direct carry-free rows:', dict(counts))


def location_checks():
    counts = Counter()
    for p in (3, 5, 7, 11):
        for d, e, j in product(range(1, 28), range(1, 28), range(1, 16)):
            g = gcd(d, e)
            if (d//g)%2 == 0 or (e//g)%2 == 0:
                continue
            s = gcd(p**d+1, p**j*(p**e+1)-1)
            assert s <= 2**d and (g == 1 or s*s <= 3**d)
            counts['resultant bounds'] += 1
            counts['nontrivial gcds'] += s > 1
    assert min(counts.values()) > 0
    print('Two unshifted pairs:', dict(counts))


if __name__ == '__main__':
    bounds()
    common_trinomials()
    four_term_checks()
    actual_row_checks()
    location_checks()
