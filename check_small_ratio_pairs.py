"""Exact common-prime check for two height-six part ratios, at bounded weights."""

import argparse
from collections import Counter, defaultdict
from itertools import combinations, combinations_with_replacement
from math import gcd, lcm

from flint import fmpz

from carry_blocks import CarryGraph, intersection
from check_part_ratio_divisors import sharp_constant
from check_primitive_row_scales import powers
from search_endpoint import avoids


def survivor_pairs():
    ratios=[pair for pair in combinations_with_replacement(range(1,7),2) if gcd(*pair)==1]
    grouped=defaultdict(list)
    for left,right in combinations_with_replacement(ratios,2):
        cap=lcm(sharp_constant(*left),sharp_constant(*right))
        if cap<=sum(left)*sum(right):
            continue
        limit=lcm(*left,*right)
        assert cap==limit and limit<2*sum(left)*sum(right)
        required=[[],[]]
        for p,w in powers(limit):
            slots=[(idx,pair) for idx,pair in enumerate((left,right)) if any(z%w==0 for z in pair)]
            assert len(slots)==1
            idx,pair=slots[0]
            unit,= [v for v in pair if v%p]
            required[idx].append((p,w,pow(unit,-1,w)))
        if len(required[0])>len(required[1]):
            left,right=right,left
            required.reverse()
        assert right in ((2,5),(3,5),(4,5),(5,6)) and len(required[1])>=2
        grouped[right].append((left,limit,*required))
    return grouped


def mixed_shapes():
    results=set()
    ratios=[pair for pair in combinations_with_replacement(range(1,7),2) if gcd(*pair)==1]
    for a,b in ratios:
        for r,s in combinations_with_replacement(range(1,7),2):
            for t in range(1,7):
                if gcd(r,s,t)>1:
                    continue
                cap=lcm(a,b,r,s,t)
                for d in range(1,cap+1):
                    if cap%d:
                        continue
                    for u in range(3,(d-1)//(a+b)+1):
                        if gcd(d,u)>1:
                            continue
                        for x in range(1,u-1):
                            for y in range(1,u-x):
                                z=u-x-y
                                if r*x+s*y!=t*z:
                                    continue
                                for i,j in combinations((x,y,z),2):
                                    h=gcd(i,j)
                                    if lcm(a,b,i//h,j//h)<=(a+b)*((i+j)//h):
                                        break
                                else:
                                    results.add((a,b,d,u,*sorted((x,y,z))))
    assert results=={(3,5,60,7,1,2,4)}
    return results


def unitary(n):
    values=[1]
    for p,e in fmpz(n).factor():
        power=int(p)**int(e)
        values+= [x*power for x in values]
    return values


def two_bases(right):
    return (3,5) if right in ((3,5),(5,6)) else (2,5)


def required_suffix(right):
    modulus,residue=1,0
    for p in two_bases(right):
        divisible,= [x for x in right if x%p==0]
        unit,= [x for x in right if x%p]
        q=p
        while divisible%(q*p)==0:
            q*=p
        residue += modulus*((pow(unit,-1,q)-residue)*pow(modulus,-1,q)%q)
        modulus *= q
    return modulus,residue


def dimension_bounds():
    upper={(2,5):{2:631,5:728},(3,5):{3:644,5:710},
           (4,5):{2:574,5:699},(5,6):{2:608,3:644,5:693}}
    length=4096
    for ratio,entries in upper.items():
        for p,numerator in entries.items():
            ways=CarryGraph(p,ratio,(0,0)).block_bound(length)
            assert ways**1000 < p**(length*numerator)
    print("digit-block exponents verified at length",length,flush=True)


def run(limit):
    dimension_bounds()
    print("mixed shapes",mixed_shapes(),flush=True)
    for right,cases in sorted(survivor_pairs().items()):
        bases=two_bases(right)
        mod,res=required_suffix(right)
        graphs=[CarryGraph(p,right,(0,0)) for p in bases]
        high=(limit-1)//sum(right)
        candidates,_,_=intersection(graphs,high,mod,res)
        if high>=1 and 1%mod==res and all(graph.accepts_integer(1) for graph in graphs):
            candidates.append(1)
        counts=Counter(scales=len(candidates))
        for h in candidates:
            divisors=unitary(h)
            for left,cap,first,second in cases:
                if gcd(h,cap)>1 or any(h%w!=r for p,w,r in second):
                    continue
                for d in divisors:
                    if d<=sum(left):
                        continue
                    t=h//d
                    z=-pow(cap*t,-1,d)%d
                    n=1+cap*t*z
                    if not 0<sum(left)*z<d or not sum(right)*h<n<=limit:
                        continue
                    g=n//d*z
                    if gcd(g,h*cap)>1 or any(g%w!=r for p,w,r in first):
                        continue
                    counts['row_pairs']+=1
                    a=(*[i*g for i in left],n-sum(left)*g)
                    b=(*[i*h for i in right],n-sum(right)*h)
                    assert min(*a,*b)>0 and sum(a)==sum(b)==n and n*(n-1)==cap*g*h
                    if any(not avoids(a,p) for p,w,r in first) or any(not avoids(b,p) for p,w,r in second):
                        continue
                    counts['small_prime_survivors']+=1
                    primes={int(p) for v in (n,n-1) for p,ex in fmpz(v).factor()}
                    if all(avoids(a,p) or avoids(b,p) for p in primes):
                        raise ArithmeticError((n,a,b))
        print(right,dict(counts),flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--limit",type=int,default=10**18)
    args=parser.parse_args()
    run(args.limit)
