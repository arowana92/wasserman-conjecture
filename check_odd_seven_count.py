"""Exact degree and two-variable gcd controls for the odd-prime counting bound."""
from fractions import Fraction
from itertools import product
from math import gcd

from flint import fmpz_poly
from check_binomial_trinomial_gcd import doubled_offsets, line_factor

X = fmpz_poly([0, 1])


def degree_controls():
    count = sharp = 0
    for p in (2, 3, 5):
        for n, e, f in product(range(1, 13), range(13), range(13)):
            if gcd(n,e,f) != 1:
                continue
            for _, s, t in doubled_offsets(p):
                for parity in (0, 1):
                    heights = (parity*n, parity*e+s, parity*f+t)
                    if any(j % 2 for j in heights):
                        continue
                    if p == 2 and len({(0,0),(e,s),(f,t)}) < 3:
                        continue
                    common = line_factor(p,(n,e,f),tuple(j//2 for j in heights))
                    size = max(n,e,f)
                    assert 2*common.degree() <= size
                    count += 1
                    sharp += 2*common.degree() == size
    assert sharp
    print(f"{count} exact trinomial degrees, including {sharp} attaining one half")


def size_checks():
    H = 10**24
    K = 1200000
    assert Fraction(3,2)+Fraction(18,10**8)+Fraction(8,H)+Fraction(4,9) < 2
    assert 50*(K+1)*32+Fraction(H+2,32) < Fraction(H,9)
    assert 2*(7*(2*K+1)*32+7*(2*K+1)*32+32) < 3000000000


def line_controls():
    # Small exact evaluations with integral powers. Use even L for
    # this independent univariate gcd test, so the coefficients are rational.
    common = coprime = 0
    samples = ((3, 0, 2, 1, 3, 0, 0), (3, 2, 3, 0, 1, 1, 1),
               (2, 4, 1, 2, 3, 1, 0), (3, 6, 1, 2, 4, 0, 1),
               (2, 5, 3, 1, 2, 0, 1))
    for p,b,n,u,v,i,j in samples:
        for q,r in product(range(1,5),range(-7,8)):
            if gcd(q,r)!=1:
                continue
            for offset in range(-8,9):
                # Evaluate the actual affine-line gcd after clearing powers
                # of T and p, using signed powers in the singleton.
                F = 1+p**(b*n)*X**(q*n)
                Y = 1+p**(b*u+i)*X**(q*u)+p**(b*v+j)*X**(q*v)
                mon = X**max(0,-r)
                scalar = p**max(0,-offset)
                R = p**max(0,offset)*X**max(0,r)
                for first,second in ((F,scalar*mon*Y+R),(Y,scalar*mon+R*F)):
                    if first.gcd(second).degree():
                        assert abs(q*offset-b*r) <= 7*(abs(r)+n*q)+4*q
                        common += 1
                    else:
                        coprime += 1
    assert common and coprime
    print(f"{common} nonconstant and {coprime} constant affine-line gcds")


if __name__ == '__main__':
    degree_controls()
    size_checks()
    line_controls()
