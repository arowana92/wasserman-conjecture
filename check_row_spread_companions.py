"""Audit the rows passing the spread search and enumerate their possible companions."""

from argparse import ArgumentParser
from collections import Counter
from functools import lru_cache
from math import gcd

from check_three_factors import valuation
from search_endpoint import find_row
from search_even_endpoint import prime_powers


@lru_cache(maxsize=65536)
def factors(n):
    return tuple(prime_powers(n).items())


@lru_cache(maxsize=65536)
def idempotents(t):
    residues = [0]
    for _,w in factors(t):
        stride = t//w
        unit = stride*pow(stride,-1,w)
        residues += [(r+unit) % t for r in residues]
    return residues


def check_by_parts(n, u, v, t, points, counts):
    # Each normalized part is zero or U modulo every full power of T.
    # Check the smaller two individually, reconstruct the third by subtraction.
    least = [(u*r) % t or t for r in idempotents(t)]
    small = [y for start in least for y in range(start,u//3+1,t) if (v*y | n) == n]
    if not small:
        return
    counts['binary-compatible smallest part'] += 1
    middle = sorted(y for start in least for y in range(start,(u-1)//2+1,t) if (v*y | n) == n)
    for y in small:
        a = v*y
        for z in middle:
            if z < y:
                continue
            third = u-y-z
            if third < z:
                break
            row = a,v*z,v*third
            if (row[0] | row[1] | row[2]) != n:
                continue
            if any(sum(k % w for k in row) != 1 for _,w in points):
                continue
            counts['full binary rows with predecessor residues'] += 1
            if all(valuation(n,row,p)==0 for p,_ in points+factors(v)):
                raise ArithmeticError(('requires independent counterexample check',n,row))


def audit(path):
    counts = Counter()
    smallest = None
    with open(path) as records:
        for line in records:
            n,u,s,t,a,d,f = map(int,line.split())
            assert n == s*t+1 and n % u == 0
            v = n//u
            assert gcd(u,v) == gcd(s,t) == 1
            row = a,a+d,a+d+f
            assert min(row)>0 and sum(row)==n and gcd(*row)%u==0 and d+f>0
            for p,_ in factors(u)+factors(s):
                assert valuation(n,row,p)==0
            counts['passing first rows'] += 1
            smallest = min(smallest or n,n)
            points = factors(t)
            if any(2*w+1 > u for _,w in points):
                continue
            counts['passing companion size test'] += 1
            check_by_parts(n,u,v,t,points,counts)
            other = find_row(n,v,factors(v),t,points)
            if other is not None:
                raise ArithmeticError(('requires independent counterexample check',n,row,other))
    print(dict(counts),'smallest weight',smallest)


def control(path):
    from search_endpoint import atoms, smallest_factors, splits

    lower, upper, spread = 35000, 55000, 1000
    spf = smallest_factors(upper)
    found = set()
    for n in range(lower,upper+1):
        predecessor = atoms(n-1,spf)
        for u,uf,_ in splits(atoms(n,spf)):
            if not 7<=u<=spread or u%2==0:
                continue
            v,h = n//u,spread//u
            for s,sf,_ in splits(predecessor):
                t = (n-1)//s
                if s%2==0 or 27*t*t > u**3:
                    continue
                for x in range(max(1,(v-2*h+2)//3),v//3+1):
                    for y in range(max(x,v-2*x-h),(v-x)//2+1):
                        row = u*x,u*y,u*(v-x-y)
                        if row[0]==row[2]:
                            continue
                        if all(valuation(n,row,p)==0 for p,_ in uf+sf):
                            found.add((n,u,s,t,row[0],row[1]-row[0],row[2]-row[1]))
    with open(path) as records:
        emitted = {r for r in (tuple(map(int,line.split())) for line in records) if lower<=r[0]<=upper}
    assert found == emitted and found
    print('Direct weight and row scan matched',len(found),'positive controls')


if __name__ == '__main__':
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('rows')
    parser.add_argument('--control',action='store_true',help='compare the 1000-spread run against direct weight enumeration')
    args=parser.parse_args()
    if args.control:
        control(args.rows)
    else:
        audit(args.rows)
