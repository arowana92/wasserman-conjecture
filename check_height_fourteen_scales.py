"""Two factorial-floor searches for the height-fourteen ratio families."""

import argparse
from bisect import bisect_right
from functools import lru_cache
from itertools import combinations, permutations, product
from math import gcd, log
from pathlib import Path

from carry_blocks import CarryGraph
from check_affine_endpoint import factors


def options(case):
    a,b,d,u,j,*row = case
    grouped = {}
    for p in sorted({2,*dict(factors(d*u))}):
        q = dict(factors(d)).get(p,1)
        left = u%p==0 or a%q==0 or b%q==0
        right = u%p!=0 and sum(v%q==0 for v in row)>=2
        grouped[p] = [i for i,yes in enumerate((left,right)) if yes]
    return grouped


def rows_at(case,n):
    a,b,d,u,j,*xs = case
    m,rem = divmod(n-1,d)
    assert not rem and n%u==0
    g = u*j*m
    return ((a*g,b*g,n-(a+b)*g),tuple(v*(n//u) for v in xs))


def prime_levels(bases,bound):
    result = {}
    for p in bases:
        seq = [p]
        while seq[-1]<bound:
            seq.append(p*seq[-1])
        result[p] = seq
    return result


def forward(case,assigned,bound):
    a,b,d,u,j,*_ = case
    modulus,residue = 1,0
    for p,w in factors(d):
        if assigned[p]==0:
            h = a if a%p else b
            residue += modulus*((pow(h*u*j,-1,w)-residue)*pow(modulus,-1,w)%w)
            modulus *= w
    m = residue+modulus*((-pow(d,-1,u)-residue)*pow(modulus,-1,u)%u)
    n = 1+d*m
    step = d*u*modulus
    start = rows_at(case,n)
    slopes = [tuple(y-x for x,y in zip(X,Y)) for X,Y in zip(start,rows_at(case,n+step))]
    levels = prime_levels(assigned,bound)
    checked = 0
    hits = []
    while n<=bound:
        parts = rows_at(case,n)
        last = {p:bisect_right(s,n) for p,s in levels.items()}
        skip = 0
        for depth in range(1,max(last.values())+1):
            for p,i in assigned.items():
                if depth>last[p]:
                    continue
                q = levels[p][last[p]-depth]
                floor = [x//q for x in parts[i]]
                if n//q>sum(floor):
                    shift = min(((k+1)*q-x+s-1)//s for x,k,s in zip(parts[i],floor,slopes[i]))
                    skip = max(skip,shift)
            if skip:
                break
        if not skip:
            hits.append(n)
        n += step*max(1,skip)
        checked += 1
    return checked,hits


def backward(case,opts,bound):
    _,_,d,u,*_ = case
    start = 1+d*(-pow(d,-1,u)%u)
    step = d*u
    n = bound-(bound-start)%step
    levels = prime_levels(opts,bound)
    checked = 0
    hits = []
    while n>=start:
        parts = rows_at(case,n)
        last = {p:bisect_right(seq,n) for p,seq in levels.items()}
        cutoffs = {p:{} for p in opts}
        low = None
        for depth in range(1,max(last.values())+1):
            for p,eligible in opts.items():
                if depth>last[p]:
                    continue
                q = levels[p][last[p]-depth]
                top = n//q
                for i in eligible:
                    if i not in cutoffs[p] and top>sum(x//q for x in parts[i]):
                        cutoffs[p][i] = top*q
                if len(cutoffs[p])==len(eligible):
                    first = max(cutoffs[p].values(),default=0)
                    low = first if low is None else min(low,first)
            if low is not None:
                break
        if low is None:
            hits.append(n)
            low = n
        n = low-1-(low-1-start)%step
        checked += 1
    return checked,hits


@lru_cache(maxsize=256)
def dimension(p,slopes,offsets):
    length = 64
    count = CarryGraph(p,slopes,offsets).block_bound(length)
    k = int(1000*log(count,p)/length)+1
    assert count**1000<p**(length*k)
    return k


def dimension_margin(case,assigned):
    a,b,d,u,j,*xs = case
    slopes = ((a*u*j,b*u*j,d-(a+b)*u*j),tuple(xs))
    dims = []
    for p,i in sorted(assigned.items(),key=lambda item:(item[0]>7,item[1]==0,item[0])):
        if p>19:
            continue
        dims.append(dimension(p,slopes[i],(0,0,int(i==0))))
        margin = 1000*(len(dims)-1)-sum(dims)
        if margin>0:
            return margin
    raise ArithmeticError(("dimension sum",case,assigned))


def minimum_plane(row):
    return min((max(r,s,t) for x,y,z in permutations(row)
                for r in range(1,15) for s in range(r,15)
                for t,rem in [divmod(r*x+s*y,z)] if not rem and 1<=t<=14),default=15)


def run(digits,rebuild):
    path = Path(__file__).with_name("height14_ratio_families.txt")
    cases = [tuple(map(int,line.split())) for line in path.read_text().splitlines()
             if line and not line.startswith("#")]
    if rebuild:
        from check_ratio_groups_cross import run as enumerate_groups
        shapes,_,survivors = enumerate_groups(14)
        assert len(shapes)==425908 and survivors==cases
    bounds = []
    totals = [0,0]
    least_margin = 1000
    assignments = 0
    for case in cases:
        *_,x,y,z = case
        row = (x,y,z)
        h = gcd(*row)
        bounds.append((case[1],minimum_plane(row),max(gcd(i,j)//h for i,j in combinations(row,2))))
        opts = options(case)
        for choices in product(*opts.values()):
            assigned = dict(zip(opts,choices))
            least_margin = min(least_margin,dimension_margin(case,assigned))
            count,hits = forward(case,assigned,10**digits)
            assert not hits
            totals[0] += count
            assignments += 1
        count,hits = backward(case,opts,10**digits)
        assert not hits
        totals[1] += count
    assert len(cases)==40 and least_margin==185
    print("families",len(cases),"assignments",assignments,"minimum dimension margin",least_margin,
          "decimal exponent",digits,"forward/backward intervals",totals,flush=True)
    for h in (6,8,10,12,14):
        planes = min([14]+[k-1 for b,k,g in bounds if b<=h])
        groups = min([14]+[g-1 for b,k,g in bounds if b<=h])
        print("ratio height",h,"certified plane height",planes,"certified pair gcd bound",groups)


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--digits",type=int,default=1000)
    parser.add_argument("--rebuild",action="store_true")
    args=parser.parse_args()
    run(args.digits,args.rebuild)
