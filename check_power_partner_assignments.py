"""Rebuild the cubic initial intervals with subset sums and factorial floors."""

from functools import lru_cache
from itertools import combinations
from math import gcd, prod

from flint import fmpz
from sympy import isprime
from check_three_factors import valuation


@lru_cache(None)
def factor(n):
    factors = tuple((int(p),int(p)**int(e)) for p,e in fmpz(n).factor())
    assert prod(w for _,w in factors) == n and all(isprime(p) for p,_ in factors)
    return factors


def prime_power(n):
    return n>=2 and len(factor(n))==1


def check_weight(n):
    assigned = rows = 0
    weight,shift = factor(n),factor(n-1)
    for length in range(3,len(shift)+1):
        for selected in combinations(shift,length):
            S = prod(w for _,w in selected)
            T = (n-1)//S
            largest = max(w for _,w in selected)
            other = max((w for _,w in shift if T%w==0),default=1)
            for width in range(1,len(weight)):
                for take in combinations(weight,width):
                    U = prod(w for _,w in take)
                    V = n//U
                    if min(U,V)<3 or V>S or V<=2*largest or U<=2*other:
                        continue
                    if 27*S*S>V**3 or 27*T*T>U**3:
                        continue
                    assigned += 1
                    # Enumerate each normalized part as a subset sum of CRT idempotents.
                    residues = [0]
                    for _,w in selected:
                        atom = S//w*(V*pow(S//w,-1,w)%w)
                        residues += [(r+atom)%S for r in residues]
                    full = len(residues)-1
                    masks = [k for k in range(1,full) if 0<residues[k]<V]
                    for i,j in combinations(masks,2):
                        if i&j:
                            continue
                        k = full^i^j
                        if k<=j or k not in masks:
                            continue
                        entries = tuple(U*residues[t] for t in (i,j,k))
                        if sum(entries)==n and all(valuation(n,entries,p)==0 for p,_ in (*take,*selected)):
                            rows += 1
    return assigned,rows


def binary_orbits():
    rejected = controls = 0
    for e in (16,24,32,48):
        P = 1<<e
        orbit = {pow(2,j,P+1) for j in range(2*e)}
        for coeffs in ((3,1,1),(1,0,2),(1,1,3),(2,1,2),(1,2,4),
                       (3,1,1,2,2),(1,2),(1,1,1),(1,2,1,1)):
            assert 8*sum(coeffs)**3 < P
            g = sum(v*P**i for i,v in enumerate(coeffs))
            if g%2==0:
                continue
            test = abs(sum(v*(-1)**i for i,v in enumerate(coeffs)))//gcd(*coeffs)
            odd_factor = test>0 and test & (test-1)!=0
            mask = (P+1)*g
            part = mask
            successes = 0
            while part:
                if part<mask and part%g==0 and part%(P+1) in orbit:
                    successes += 1
                part = (part-1)&mask
            if odd_factor:
                assert successes==0
                rejected += 1
            else:
                controls += successes
    assert rejected and controls
    print(f"{rejected} binary polynomial exclusions and {controls} admissible "
          "part controls in the other examples.")


def verify():
    box = 32
    weights = set()
    largest = 0
    # A prime-power pair containing consecutive integers contains an even power of two.
    even_powers = [1<<e for e in range(1,(1+4*16*box).bit_length()*5+2)]
    inputs = {v for t in even_powers for v in (t-1,t) if prime_power(v) and prime_power(v+1)}
    for x in inputs:
        for m in range(1,box+1):
            for h in range(1,box+1):
                z = m*x+h
                if not prime_power(z) or gcd(z,x*(x+1))!=1:
                    continue
                for d in (1,2,3,4,6,8,12,16):
                    if x>(1+2*d*(m+h))**5:
                        continue
                    n = 1+d*x*(x+1)*z
                    largest = max(largest,n)
                    if n>10**9 and len(factor(n-1))>=3:
                        weights.add(n)
    assignments = rows = 0
    for n in weights:
        total,valid = check_weight(n)
        assignments += total
        rows += valid
    assert (len(weights),assignments,rows,largest)==(3134,27480,0,5070602398182799483249597874177)
    print(f"{len(weights)} weights beyond the initial interval, {assignments} assignments, "
          f"{rows} surviving rows, largest weight {largest}.")


if __name__=="__main__":
    binary_orbits()
    verify()
