"""Classify a second-interval row when its smaller predecessor power is at most 128."""

from collections import Counter
from functools import lru_cache
from itertools import combinations
from math import gcd, isqrt, lcm

from carry_blocks import CarryGraph, intersection
from check_endpoint_pair_bound import prime_powers
from check_shifted_directions import possible_scales
from check_shifted_prefixes import floor_cover
from search_endpoint import atoms, smallest_factors


@lru_cache(maxsize=None)
def divisors(n):
    return tuple(i for i in range(1,n+1) if n%i == 0)


@lru_cache(maxsize=None)
def triples(total):
    return tuple((f,g,total-f-g) for f in range(1,total//3+1)
                 for g in range(f,(total-f)//2+1))


@lru_cache(maxsize=None)
def by_cofactor(u, scale):
    result = set()
    for r in range(scale+1, isqrt(u*u*scale//27)+1, scale):
        forced = [w for _,w in prime_powers(r) if (u//scale)%w]
        if any(2*w >= scale for w in forced):
            continue
        for other in triples(scale):
            if all(sum(e%w == 0 for e in other) >= 2 for w in forced):
                result.add((r,other))
    return result


@lru_cache(maxsize=None)
def by_companion(u, scale):
    result = set()
    for other in triples(scale):
        available = lcm(*(gcd(a,b) for a,b in combinations(other,2)))
        for t in divisors(available):
            for w in divisors(u//scale):
                r = w*t
                if r > 1 and gcd(w,t) == 1 and r%scale == 1:
                    result.add((r,other))
    return result


def classify(limit):
    # r <= (P-1)^(3/2) / sqrt(27).
    spf = smallest_factors(limit*(1+isqrt(limit**3//27)))
    counts, survivors = Counter(), []
    for P in range(2,limit+1):
        if len(tuple(prime_powers(P))) != 1:
            continue
        for u in range(3,P):
            if gcd(P,u) != 1:
                continue
            for scale in divisors(u):
                if scale < 3:
                    continue
                cases = by_cofactor(u,scale)
                assert cases == by_companion(u,scale), (u,scale)
                for r,other in sorted(cases):
                    if gcd(r,P) != 1:
                        continue
                    counts['companion choices'] += 1
                    m = r*P
                    a, b = sorted((u,(r-1)*P))
                    shape = (scale,a,b,P-u,*other)
                    for parity, (residue,step), primes in possible_scales(m,atoms(m,spf),shape,spf):
                        survivors.append((P,u,scale,r,other,parity,residue,step,tuple(primes)))
    return counts, survivors


def check():
    counts, surviving = classify(128)
    print(dict(counts))
    expected = [(101,77,7,22,(1,2,4),1,8101,15554,(7,11,101)),
                (127,104,4,13,(1,1,2),1,2297,6604,(2,13,127))]
    assert surviving == expected, surviving
    for P,u,scale,r,other,parity,residue,step,bases in surviving:
        slopes = (u,(r-1)*P,P-u)
        graphs = [CarryGraph(p,slopes,(0,0,1)) for p in bases]
        found, nodes, _ = intersection(graphs,10**1000,step,residue)
        direct, intervals = floor_cover(slopes,bases,step,residue,10**1000)
        assert direct == found == []
        print('P, inverse, scale, cofactor, normalized companion, parity, progression, bases:',
              (P,u,scale,r,other,parity,residue,step,bases))
        print('Graph prefix nodes and separate floor intervals through 10^1000:', nodes, intervals)
    orbits = []
    for height in range(2021,2122):
        fs = list(prime_powers(height))
        if len(fs) != 1:
            continue
        q = fs[0][0]
        pos = {}
        value = 1
        while value not in pos:
            pos[value] = len(pos)
            value = value*q%15554
        if 8101 in pos:
            assert value == 1
            orbits.append((q,height,pos[8101],len(pos)))
    assert orbits == [(2081,2081,214,300)]


if __name__ == '__main__':
    check()
