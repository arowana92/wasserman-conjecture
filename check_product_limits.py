"""Check the limiting family recurrence and the finite binary exclusions."""

from itertools import combinations_with_replacement
from math import comb, gcd, prod

from sympy import Poly, Rational, expand, symbols

n = symbols('n')


def verify_degrees(f, g, first_bin, deficits):
    f1,f2,f3 = f
    g2,g3 = g
    d0,d1,d2,d3 = deficits
    a = (f1*n+d0)/first_bin
    b = (f2*a*n+d1)/d0
    c = (f3*a*b*n+d2)/d1
    v = f1*b*c+f2*a*c+f3*a*b
    t = (g2*first_bin*v+d3)/d2
    w = (g3*first_bin*t*v-1)/d3
    u = n*t*w+g2*first_bin*w+g3*first_bin*t
    m = a*b*c*first_bin*t*w
    assert expand(u*v-m-1) == 0
    poly = Poly(m,n)
    denom = first_bin**9*d0**8*d1**4*d2**2*d3*f1**9*f2**17*f3**21*g2**23*g3**24
    assert poly.degree()==25 and poly.LC()/prod(f+g)**25 == Rational(1,denom)


def binary_exclusions():
    for last,bits in ((3,5),(7,9),(9,9)):
        mod = 2**bits
        for v0 in range(0,mod,2):
            a = v0+1
            b = v0*a+1
            c = v0*a*b+1
            v = a*b+a*c+b*c
            d = v+last
            e = ((v*d-1)*pow(last,-1,mod)) % mod
            row = [(v*z)%mod for z in (v0*d*e,e,d)]
            assert sum(row) > (row[0]|row[1]|row[2])
        print(f'Last deficit {last}: forced a binary carry in all {mod//2} even residues modulo {mod}.')


def squarefree_denominators():
    modulus = 512
    excluded = 0
    for last in range(1,modulus,2):
        inverse = pow(last,-1,modulus)
        for n0 in range(0,modulus,2):
            a = n0+1
            b = n0*a+1
            c = n0*a*b+1
            scale = (a*b+a*c+b*c) % modulus
            d = scale+last
            e = (scale*d-1)*inverse % modulus
            x,y,z = (scale*t % modulus for t in (n0*d*e,e,d))
            if not (x&y or x&z or y&z):
                break
        else:
            excluded += 1
    assert excluded == 99
    print(f'Excluded {excluded} of 256 odd last-deficit classes modulo {modulus}.')


def exact_two_rows(bound):
    total = odd_common = 0
    for weight in range(6,bound+1):
        rows = []
        predecessor = weight-1
        for a in range(2, weight//3+1):
            for b in range(a,(weight-a)//2+1):
                c = weight-a-b
                coefficient = comb(weight,a)*comb(weight-a,b)
                rows.append((a*b*c, gcd(weight,coefficient), gcd(predecessor,coefficient)))
        for (a,d,j),(b,e,h) in combinations_with_replacement(rows,2):
            d0 = gcd(d,e)
            j0 = gcd(j,h)
            value, remainder = divmod(d0**3*j0**2*a*b,weight**3*predecessor**2)
            assert remainder == 0 and value >= 1
            common = d0*j0
            assert (predecessor//j0)*common**9 < 9**40*value**25
            total += 1
            if common % 2:
                odd_common += 1
                assert value % (8 if weight % 2 == 0 else 2) == 0
    print(f'Checked the unrestricted inequality on {total:,} exact pairs, including {odd_common:,} with odd endpoint gcd.')


if __name__ == '__main__':
    exact_two_rows(55)
    for f,g,t,d in [
        ((1,1,1),(1,1),1,(1,1,1,1)),
        ((2,1,3),(1,2),5,(1,4,2,7)),
        ((1,4,2),(3,5),3,(2,1,5,4)),
        ((2,3,5),(7,11),2,(3,5,7,11)),
    ]:
        verify_degrees(f,g,t,d)
    binary_exclusions()
    squarefree_denominators()
