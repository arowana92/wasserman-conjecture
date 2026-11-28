"""Compute the least E with uV-1 dividing E times the coefficient of Vy."""

import argparse
from functools import lru_cache
from heapq import heappop, heappush
from itertools import combinations
from math import gcd, lcm

from check_affine_endpoint import factors
from check_direction_pair_count import heavy_rows, valuation


def vp(n,p):
    e=0
    while n%p==0:
        e+=1
        n//=p
    return e


def decrement_witness(row,p,e,j):
    """Find an actual multiplier attaining defect at least j, if one exists."""
    q=p**j
    r=pow(sum(row),-1,q)
    unit,= [i for i,x in enumerate(row) if x%p]
    vanishing=next(i for i,x in enumerate(row) if x%p==0 and vp(x,p)==e)
    terms=tuple(x*r-int(i in (unit,vanishing)) for i,x in enumerate(row))
    assert sum(x%q for x in terms)==q-1
    start=(*[x//q for x in terms],0)
    finish=(0,0,0,0)
    budget=e-j
    dist={start:0}
    queue=[(0,start)]
    previous={}
    while queue:
        cost,state=heappop(queue)
        if cost!=dist[state]:
            continue
        if state==finish:
            n=0
            while state!=start:
                state,digit=previous[state]
                n=p*n+digit
            return r+q*n
        for digit in range(p):
            values=[offset+s*digit for s,offset in zip(row,state[:3])]
            carry=(state[3]+sum(v%p for v in values))//p
            dest=(*[v//p for v in values],carry)
            new=cost+carry
            if new<=budget and new<dist.get(dest,budget+1):
                dist[dest]=new
                previous[dest]=state,digit
                heappush(queue,(new,dest))
    return None


@lru_cache(maxsize=100000)
def prime_loss(row,p,w):
    e=vp(w,p)
    for j in range(e,0,-1):
        v=decrement_witness(row,p,e,j)
        if v is not None:
            n=sum(row)*v
            actual=vp(n-1,p)-valuation(n,tuple(v*x for x in row),p)
            assert actual==j
            return p**j,v
    return 1,None


def constant(row):
    if min(row)<1 or gcd(*row)!=1:
        raise ValueError("a positive primitive triple is required")
    coarse=lcm(*(gcd(a,b) for a,b in combinations(row,2)))
    answer=1
    for p,w in factors(coarse):
        power,_=prime_loss(row,p,w)
        answer*=power
    return answer


def census(limit):
    # H must reach u before the exact constant can do so.
    possible=sorted({row for row,_ in heavy_rows(limit)},key=lambda y:(sum(y),y))
    survivors=[]
    for row in possible:
        c=constant(row)
        if c>=sum(row):
            survivors.append((row,c))
    print("sum limit",limit,"coarse heavy directions",len(possible),
          "exact heavy directions",len(survivors),"first",survivors[:5],flush=True)
    return survivors


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--row",type=int,nargs=3)
    parser.add_argument("--limit",type=int,default=256)
    args=parser.parse_args()
    if args.row:
        row=tuple(args.row)
        print("sharp constant",constant(row))
    else:
        census(args.limit)
