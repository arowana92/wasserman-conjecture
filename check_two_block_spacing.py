"""Check the divisibilities behind exponent spacing, including both proportionality cases."""

from collections import Counter
from itertools import combinations, product
from math import gcd

from check_three_factors import powers, valuation


def splits(n, p):
    return [(x,y,n-x-y) for x in range(n+1) for y in range(n-x+1)
            if valuation(n,(x,y,n-x-y),p) == 0]


def check_shape(p, q, b, c, upper, lower, atoms, counts):
    d = b*q+c
    parts = tuple(i*q+j for i,j in zip(upper,lower))
    if 0 in parts:
        return
    assert valuation(d,parts,p) == 0
    det = [b*j-c*i for i,j in zip(upper,lower)]
    if det[0] == det[1] == 0:
        counts['fully proportional'] += 1
        assert det[2] == 0 and sum(upper) == b and min(upper) > 0
        g = gcd(b,c)
        assert all(x%(d//g) == 0 for x in parts)
        return
    i = 0 if det[0] else 1
    x,y = parts[i],parts[1-i]
    mixed = det[1-i] == 0
    records = []
    for a in range(1,11):
        t = p**a
        n = d*t+1
        actual = (t*x,t*y,t*parts[2]+1)
        u = gcd(n,gcd(x,y))
        s = f = 1
        for w in atoms:
            residues = [z%w for z in actual]
            if sorted(residues) == [0,0,1]:
                s *= w
                if residues[0] == 1:
                    f *= w
        assert gcd(d,x*y)%s == 0
        if mixed:
            assert upper[1-i]%u == 0
            assert x%(s//f) == 0 and s//f <= b*c
        else:
            assert abs(det[0]*det[1])%s == 0
        records.append((a,u,f))
    for (a,u,f),(ap,up,fp) in combinations(records,2):
        left, right, ambient = (f,fp,d) if mixed else (u,up,x)
        shared = gcd(left,right)
        assert ambient*shared >= left*right
        assert pow(p,ap-a,shared) == 1%shared
        if shared > 1:
            counts['shared predecessor factor' if mixed else 'shared weight factor'] += 1


def main():
    counts = Counter()
    for p,b,c in product((2,3,5),range(1,9),range(1,8)):
        q = p
        while q <= c:
            q *= p
        for q0 in (q,p*q):
            atoms = tuple(powers(b*q0+c).values())
            for upper,lower in product(splits(b,p),splits(c,p)):
                check_shape(p,q0,b,c,upper,lower,atoms,counts)
    assert len(counts) == 3 and min(counts.values()) > 0
    print(dict(counts))


if __name__ == '__main__':
    main()
