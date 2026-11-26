"""Exact finite comparisons for the seven-unit large-span theorem at 5 and 7."""
from collections import Counter
from fractions import Fraction
from itertools import combinations
from math import gcd

from flint import fmpz_poly, nmod_poly
from sympy import Poly, Symbol

from check_gap_resultants import resultant
from check_seven_refined_cutoffs import below_one

X = fmpz_poly([0, 1])


def constants():
    minimum = 10**16
    assert 200000**3 <= minimum and 10000**4 == minimum
    assert Fraction(18, 200000)+Fraction(14, 10000)+Fraction(8, minimum) < Fraction(1, 500)
    assert Fraction(3, 2)+Fraction(1, 500)+Fraction(26, 200000) < 2
    assert Fraction(14, 15)+Fraction(1, 500)+Fraction(1700, 10**8) < 2
    for p, main, other in (
        (5, Fraction(3, 5), [(3, 2, 1)]),
        (5, 0, [(4, 2, 1)]),
        (5, Fraction(22, 35), [(3, 2, 1)]),
        (7, Fraction(6, 7), [(3, 2, 1)]),
        (5, Fraction(14, 15), [(2, 20, 9)]),
        (5, Fraction(7, 10), [(8, 1, 1)]),
    ):
        power = Fraction(main) + Fraction(1, 500) - 2
        assert below_one([(p, power.numerator, power.denominator), *other])
    print("All uniform logarithmic comparisons hold by rational arithmetic")


def repeated_resultants():
    checks = 0
    for w in range(1, 16):
        for n, ell in ((3, 1), (1, 3), (2, 3), (3, 3)):
            for c in range(-3*w, 5*w+1):
                h = max(3*w, max(c, ell*w)-min(0, c))
                g = gcd(w, c)
                low = min(0, c)
                f = 2+X**(n*w//g)
                other = sum(X**((a-low)//g) for a in (0, ell*w, c))
                if f == other:
                    continue
                r = resultant(f, other)
                assert 0 < r and r**(9*g) <= 8**(9*g)*2**(9*h+3*w)
                checks += 1
    print(f"{checks} exact repeated-position resultant bounds")


def gram_controls():
    count = 0
    shapes = [1+X**j+X**d for d in range(2, 23) for j in range(1, d)]
    for f, h in combinations(shapes, 2):
        height = max(f.degree(), h.degree())
        r = resultant(f, h)
        assert r*r <= 8**height
        count += 1
    print(f"{count} distinct-position determinant bounds")


def normalized_profiles():
    # Enumerate every four-position subset after translation and scaling,
    # independent of the matching construction in the proof.
    counts = Counter()
    for n in range(5, 40, 2):
        if n % 3 == 0:
            continue
        degree = 3*n
        limit = (20*(n+2)-1)//7
        binomial = 1+X**degree
        for a, b, c in combinations(range(1, limit+1), 3):
            places = (0, a, b, c)
            # Each opposite matching must have one pair through zero.
            # Test the two possible indices separately and retain all
            # three choices of partner. Confirm the actual gcd below.
            if not any(j % (2*n) == n and (k-l) % (2*n) == n
                       for j, k, l in ((a,b,c), (b,a,c), (c,a,b))):
                continue
            if not any(j % 6 == 3 and (k-l) % 6 == 3
                       for j, k, l in ((a,b,c), (b,a,c), (c,a,b))):
                continue
            quad = sum(X**i for i in places)
            if binomial.gcd(quad).degree() != n+2:
                continue
            counts['profiles at 5'] += 1
            at_seven = n == 5 and c <= 16
            counts['profiles at 7'] += at_seven
            for single in places:
                other = [i for i in places if i != single]
                y = sum(X**(i-min(other)) for i in other)
                mod5 = lambda f: nmod_poly([int(t) for t in f.coeffs()], 5)
                assert mod5(binomial).gcd(mod5(y)).degree() == 0
                counts['base 5 modular gcds'] += 1
                if at_seven:
                    x = Symbol("x")
                    exact = [Poly.from_list(list(reversed(f.coeffs())), x) for f in (binomial, y)]
                    assert resultant(binomial, y) == 271 == abs(exact[0].resultant(exact[1]))
                    counts['base 7 resultants'] += 1
    assert counts['profiles at 7'] == 2 and counts['base 7 resultants'] == 8
    # A completely unfiltered pass over the small base-7 case.
    admissible = []
    for rest in combinations(range(1, 17), 3):
        places = (0, *rest)
        if (1+X**15).gcd(sum(X**i for i in places)).degree() == 7:
            admissible.append(places)
    assert len(admissible) == 2
    print(dict(counts), "base 7 positions", admissible)


if __name__ == '__main__':
    constants()
    repeated_resultants()
    gram_controls()
    normalized_profiles()
