"""Verify the 31 and mixed 11 cutoffs for seven freely adding digit units."""
from fractions import Fraction
from itertools import combinations, combinations_with_replacement
from math import gcd

from check_gap_resultants import resultant
from check_seven_large_bases import X, poly
from check_seven_refined_cutoffs import below_one


def powers():
    b = 31
    # Each row bounds (U^3 A^2 C^2)^(1/H), with H >= 101.
    rows = [
        [(3,3,2),(6,1,1),(9,1,1)],                                     # k = 0, distinct
        [(32,2,3),(b+1,2,3),(9,1,1)],                                   # L <= 2
        [(3,3,2),(6,1,1),(b+1,2,101),(9,1,1)],                          # L > 2, k = 1
        [(3,3,4),(6,1,2),(b+1,2,3),(9,1,1)],                            # divide positions by k
        [(Fraction(19,4),2,5),(4,4,15),(b,14,15),(9,1,1)],               # two matches, g >= 2
        [(Fraction(133,33),99,100),(4,33,50),(b,17,25),(9,1,1)],         # g = 1, d >= 300
        [(Fraction(27,8),1,101),(8,1,1),(10,1,1),(9,1,1)],             # two residues, k = 0
        [(Fraction(27,8),1,101),(8,1,1),(6,1,1),(b+1,2,101),(9,1,1)],
        [(5,3,4),(6,1,2),(b+1,2,3),(9,1,1)],                            # two residues, k >= 2
        [(Fraction(27,8),1,101),(8,1,1),(16,1,1),(6,1,1),(b+1,2,101)], # C > 3^H
        [(3,3,2),(54,1,1),(b,12,101)],                                  # nontrivial h_U, t <= 2
        [(3,7,18),(54,1,1),(b,2,3)],                                    # t >= 3
    ]
    for i, row in enumerate(rows):
        assert below_one(row+[(b,-2,1)]), i
    b = 11
    mixed = [
        [(3,1,1),(10,1,2)],
        [(4,1,1),(6,1,2),(b+1,1,101)],
        [(3,1,1),(2,2,3),(b+1,1,3)],
        [(3,1,1),(6,1,4),(b+1,1,3)],
        [(3,1,1),(2,8,45),(b,7,15)],
        [(3,1,1),(2,10,33),(b,13,33)],
        [(3,1,1),(2,21,101),(b,98,505)],
        [(3,1,1),(2,33,50),(b,17,50)],
    ]
    for i, row in enumerate(mixed):
        assert below_one(row+[(2,1,101),(b,-1,1)]), i
    print(f"{len(rows)} general and {len(mixed)} mixed rational-power inequalities")


def small_matching_indices():
    remaining = []
    for n in range(5,100,2):
        if n % 3 == 0:
            continue
        d, k = 3*n, n+2
        h = max(101,d)
        # The ratios decrease when H increases, and when the base increases.
        assert (Fraction(3*d-k,d-k)**(3*(d-k))*4**(2*(d-k))*31**(4*k)*9**(2*h)
                < 729*31**(4*h))
        for b in (11,12):
            h = max(101,d)
            while 2*3**h * 2**(d-k)*b**k >= b**h:
                remaining.append((b,d,h))
                h += 1
    assert remaining == [(11,105,105)]
    f = 1+11**105
    count = 0
    for a, b in combinations(range(71),2):
        pts = (a,a+35,b,b+35)
        if len({v % 105 for v in pts}) != 4:
            continue
        total = sum(11**i for i in pts)
        acap = gcd(f,total)
        for i in pts:
            z = 11**i
            ccap = gcd(total-z,f+z)
            assert 2*acap*ccap < f+total
            count += 1
    print(f"{count} direct allocations for the residual mixed case")


def companion_resultants():
    count = 0
    for d in range(1,18):
        for e, f in combinations(range(1,22), 2):
            for a in (X**d+2, 2*X**d+1):
                b = 1+X**e+X**f
                r = resultant(a,b)
                assert 0 < r <= 3**max(d,f)
                count += 1
    print(f"{count} distinct/repeated trinomial resultant comparisons")


def matching_controls():
    checked = 0
    for d in (15,21,33,35,45):
        f = 1+X**d
        # Include negative signed residues by taking representatives modulo 2d.
        for extra in combinations_with_replacement(range(2*d), 3):
            a = (0, *extra)
            matches = set()
            for j in range(1,4):
                l, r = (i for i in range(1,4) if i != j)
                u, v = a[j], abs(a[l]-a[r])
                g = gcd(d,u,v)
                if (d//g)%2 and (u//g)%2 and (v//g)%2:
                    matches.add(g)
            if not matches or d in matches:
                continue
            if len(matches) == 1:
                k = next(iter(matches))
            else:
                r,s = sorted(matches)
                assert len(matches) == 2
                k = r+s-gcd(r,s)
            if 3*k <= d:
                continue
            assert len(matches) == 2
            g = gcd(*matches)
            assert d%(3*g) == 0 and k == d//3+2*g
            assert all(i%g == 0 for i in a)
            assert f.gcd(poly(a)).degree() == k
            checked += 1
    print(f"{checked} four-unit profiles with common degree above one third")


if __name__ == '__main__':
    powers()
    small_matching_indices()
    companion_resultants()
    matching_controls()
