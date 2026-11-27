"""Exact bounds and local isometry controls for moving prime-power parameters."""
from fractions import Fraction
from math import gcd

from check_cofactor_986_powers import ROWS, productive_edges
from check_direction_prime_bases import odd_shapes, exponent_growth


def valuation(n,p):
    assert n > 0
    v = 0
    while n%p == 0:
        n //= p
        v += 1
    return v


def path_bound(slopes,numerator):
    edges = productive_edges(7,slopes)
    weights = {s:1 for s in edges}
    for _ in range(200):
        weights = {s:weights[s]+sum(weights[t] for t in row.values()) for s,row in edges.items()}
    assert all(sum(weights[t] for t in row.values())**1000 <= 7**numerator*weights[s]**1000
               for s,row in edges.items())


def shifted_digits(p,q,step,scale,offset,denom,depth,geometric):
    """Count a full block modulo p**depth after a possibly very long fixed prefix."""
    width = p**depth
    assert gcd(q,p) == gcd(denom,p) == 1
    if geometric:
        # Unit ratio congruent to one mod p (modulo 4 for p=2).
        start = 17
        stride = step*(2 if p==2 else p-1)
        if p == 2:
            h = valuation(q*q-1,2)+valuation(stride//2,2)
        else:
            h = valuation(q**(p-1)-1,p)+valuation(stride//(p-1),p)
        h += valuation(scale,p)
        modulus = p**(h+depth)
        values = [pow(q,start+stride*i,modulus) for i in range(width)]
    else:
        exponent = p**32*17+2*p**32
        stride = step*(4 if p==2 else p)
        h = valuation(stride,p)+valuation(exponent,p)+valuation(scale,p)
        modulus = p**(h+depth)
        values = [pow(q+stride*i,exponent,modulus) for i in range(width)]
    values = [((v*scale+offset)*pow(denom,-1,modulus))%modulus for v in values]
    common = p**h
    assert len({v%common for v in values}) == 1
    assert sorted(v//common for v in values) == list(range(width))


def verify():
    upper = Fraction(733,500)
    assert upper**3 > upper**2+1 and upper**125 < 2**69
    assert Fraction(1,2)+Fraction(9,10)*Fraction(69,125) < Fraction(997,1000)
    for (slopes,_),t in zip(ROWS,(284,248)):
        path_bound(slopes,t)
    for p,q in ((2,3),(3,5),(5,7),(7,11)):
        for geometric in (False,True):
            shifted_digits(p,q,p**16,p**21,1,11 if p!=11 else 13,3,geometric)
    cases = 0
    for d,shape,_,_,bases in odd_shapes():
        if d in (364,880):
            continue
        exponent_growth(d,shape,bases)
        cases += 1
    print("Uniform exponents 0.336 for bounded cofactors, 0.284 and 0.248 at 986.")
    print(cases,"odd-scale rows,",8,"complete offset-digit permutation controls.")


if __name__ == '__main__':
    verify()
