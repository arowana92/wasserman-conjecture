"""Exact checks for denominator cancellation and the two-target quadratic case."""

from collections import Counter
from fractions import Fraction
from itertools import combinations, product
from math import gcd, lcm


def common_denominators():
    counts = Counter()
    for p, g, d, s in product((2,3,5), (4,7), range(2,8), range(1,7)):
        if s >= d or gcd(d,s) != 1:
            continue
        for t in (Fraction(1), Fraction(3,5), Fraction(5,7), Fraction(2,3)):
            if t.numerator % p == 0 or t.denominator % p == 0:
                continue
            for j,k,sign in product((-1,0,2),(-2,0,1),(-1,1)):
                a,b = j+g*d,k+g*s
                x = sign**d * Fraction(p)**j * t**d
                y = sign**s * Fraction(p)**k * t**s
                den = lcm(x.denominator,y.denominator)
                i, z = int(den*x), int(den*y)
                assert gcd(den,i,z) == 1
                for h in (1,3,5,15):
                    delta, left, right = h*den, h*i, h*z
                    w = gcd(delta*p**a-left,delta*p**b-right)
                    while w % p == 0:
                        w //= p
                    wp = w//gcd(w,h)
                    assert gcd(den,wp) == gcd(t.denominator*t.numerator,wp) == 1
                    assert (t.denominator*p**g - sign*t.numerator) % wp == 0 or (
                        (t.denominator*p**g + sign*t.numerator) % wp == 0)
                    assert w <= delta * (t.denominator*p**g + t.numerator)
                    counts["modular reconstructions"] += 1
                    if gcd(delta,w)>1 and wp>2:
                        counts["noninvertible original denominator"] += 1
    return counts


def two_targets():
    coords = [v for v in product(range(5),range(9),range(6)) if v[0] or v[1]]
    targets = [(0,0,0),(1,0,0),(0,1,0),(0,0,1)]
    hits = 0
    for q in (Fraction(1,2),Fraction(1),Fraction(2)):
        for u,v in combinations(map(Fraction,("1/2","1","3/2","2","3")),2):
            choices = {(i,j):[] for i,j in product(range(2),repeat=2)}
            for d,m,c in coords:
                pair = (d*q*u*u-m*u+c, d*q*v*v-m*v+c)
                if pair in choices:
                    choices[pair].append((d,m,c))
            for first,second in product(targets,repeat=2):
                if first == second:
                    continue
                for matrix in product(*(choices[i,j] for i,j in zip(first,second))):
                    D,M,C = zip(*matrix)
                    assert sum(x>0 for x in D) >= 2
                    assert all(M)
                    assert sum(x>0 for x in C) >= 2
                    hits += 1
    assert hits
    for b in range(7,151):
        t = 2**b
        parts = [d*t*t+m*t+c for d,m,c in ((1,2,1),(4,12,8),(32,96,64))]
        assert sum(parts) == 37*t*t+110*t+73
        assert not any(i & j for i,j in combinations(parts,2))
        assert gcd(*parts) % (t+1) == gcd(parts[0]-1,*parts[1:]) % (t+2) == 0
    return hits


if __name__ == "__main__":
    assert (2-Fraction(5,3)-Fraction(9,64))*16384 > 1185+14
    print("Modular checks:", dict(common_denominators()))
    print("Two-target rational examples:", two_targets())
