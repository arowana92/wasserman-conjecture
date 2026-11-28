"""Compute the least C with N-1 dividing g*C*K on every row (ag,bg,z)."""

import argparse
from functools import lru_cache
from heapq import heappop, heappush
from itertools import combinations_with_replacement
from math import gcd

from check_line_divisors import valuation
from check_primitive_row_scales import powers


def minimum_decrement(a, b, p, full):
    """Minimize the binomial valuation on (a*g-1,b*g-1) over p-unit g."""
    unit = a if a%p else b
    residue = pow(unit,-1,full)
    for g in range(residue,residue+12*full,full):
        if valuation((a*g-1,b*g-1),p)==0:
            return 0,g
    start = (-1,-1,0)
    end = (0,0,0)
    parent = {}
    best = {start:0}
    queue = [(0,start)]
    while queue:
        cost,state = heappop(queue)
        if best[state]!=cost:
            continue
        if state==end:
            g=0
            while state!=start:
                state,digit=parent[state]
                g = p*g+digit
            assert g%p and valuation((a*g-1,b*g-1),p)==cost
            return cost,g
        x,y,carry=state
        for digit in range(int(state==start),p):
            u,i=divmod(x+a*digit,p)
            v,j=divmod(y+b*digit,p)
            more=(i+j+carry)//p
            dest=(u,v,more)
            if dest not in best or cost+more<best[dest]:
                best[dest]=cost+more
                parent[dest]=state,digit
                heappush(queue,(cost+more,dest))
    raise AssertionError("any unit multiplier can be completed")


@lru_cache(maxsize=None)
def sharp_constant(a,b):
    if min(a,b)<1 or gcd(a,b)!=1:
        raise ValueError("positive coprime terms required")
    result=1
    for p,full in powers(a*b):
        missing,_ = minimum_decrement(a,b,p,full)
        result *= max(1,full//p**missing)
    return result


def census(height):
    pairs=[(a,b) for a,b in combinations_with_replacement(range(1,height+1),2)
           if gcd(a,b)==1]
    reduced=[]
    strong=[]
    for a,b in pairs:
        c=sharp_constant(a,b)
        if c<a*b:
            reduced.append((a,b,c))
            if c<=a+b and a>1:
                strong.append((a,b,c))
    print(f"height={height}, ratios={len(pairs)}, reduced={len(reduced)}, noninteger ratios with C<=sum={len(strong)}")
    print("noninteger ratios with C<=sum:",strong)
    return reduced


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ratio",type=int,nargs=2,metavar=("A","B"))
    parser.add_argument("--height",type=int,default=128)
    args=parser.parse_args()
    if args.ratio:
        a,b=args.ratio
        print(f"C={sharp_constant(a,b)}")
        for p,full in powers(a*b):
            m,g=minimum_decrement(a,b,p,full)
            print(f"p={p}, minimum={m}, witness multiplier={g}")
    else:
        census(args.height)
