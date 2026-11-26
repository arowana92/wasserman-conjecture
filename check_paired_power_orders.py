"""Check the paired carry lemmas on exact binomial coefficients and verify both digit bounds."""

from fractions import Fraction
from math import comb, gcd

from flint import fmpz_mat

from carry_blocks import CarryGraph


def digit_matrix(base, offset, first=3, second=2):
    # Carry pair for the expressions first*x and second*x+offset.
    states = [(0,offset)]
    edges = []
    for u,v in states:
        out = []
        for digit in range(base):
            a,b = first*digit+u, second*digit+v
            if a%base + b%base < base:
                child = (a//base,b//base)
                out.append(child)
                if child not in states:
                    states.append(child)
        edges.append(out)
    return fmpz_mat([[out.count(child) for child in states] for out in edges])


def check_counts():
    # Constant coefficients come first in FLINT.
    binary = digit_matrix(2,0)
    quinary = digit_matrix(5,1)
    assert list(binary.charpoly()) == [-1,1,-1,-1,1]
    assert list(quinary.charpoly()) == [0,-1,2,1,-3,1]
    x,y = Fraction(1513,1000), Fraction(2247,1000)
    assert x**4-x**3-x**2+x-1 > 0 and y**3-2*y**2-y+1 > 0
    assert x**1000 < 2**598 and y**1000 < 5**504
    print('Verified the binary and quinary matrices and rational exponent bounds.')


def check_three_base_dimensions():
    length = 256
    for w,a,bases,expected in ((5,3,(2,3,5),(61,69,51)),
                               (7,3,(2,3,7),(57,66,62)),
                               (7,5,(2,5,7),(65,74,59))):
        bounds = []
        for p in bases:
            count = CarryGraph(p,(a,w-a),(0,1)).block_bound(length)
            powered = digit_matrix(p,1,a,w-a)**length
            assert count == sum(powered.entries())
            bounds.append(next(n for n in range(101) if count**100 <= p**(length*n)))
        assert tuple(bounds)==expected and sum(bounds)<200
        print('Factor, multiplier, bases, dimension hundredths:',w,a,bases,bounds)


def check_binomial_layers():
    cases = orders = gaps = carry_free = 0
    for p in (2,3,5,7,11):
        for h in range(1,10):
            power = p**h
            for d in range(2,201):
                if gcd(p,d)>1:
                    continue
                q = power//d+1
                if q>150 or q%p==0:
                    continue
                c = d*q-power
                residue = 1
                near = []
                returning = False
                for j in range(h):
                    near.append((d-residue,p**(h-j)))
                    residue = residue*p%d
                    returning |= residue==1
                for w in range(2,min(15,d//c)+1):
                    for a in range(1,w//2+1):
                        cases += 1
                        gap = any(w*t<=d and c<t*m for t,m in near)
                        divisible = comb(w*q,a*q)%p==0
                        if returning or gap:
                            assert divisible
                        orders += returning
                        gaps += gap
                        if not divisible:
                            carry_free += 1
                            radix = p
                            while radix < w:
                                radix *= p
                            x = q
                            while x:
                                x,digit = divmod(x,radix)
                                assert digit < radix-radix//w
    assert min(orders,gaps,carry_free)>0
    print('Binomial cases, order obstructions, gap obstructions, carry-free cases:',
          cases,orders,gaps,carry_free)


if __name__=='__main__':
    check_counts()
    check_three_base_dimensions()
    check_binomial_layers()
