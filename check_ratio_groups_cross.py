"""Enumerate denominator powers by CRT and test relaxed one-prime languages."""

import argparse
from collections import Counter
from functools import lru_cache
from itertools import combinations, combinations_with_replacement, permutations, product
from math import gcd, lcm, prod

from check_affine_endpoint import factors, divisors
from check_part_ratio_divisors_cross import transport_cap


@lru_cache(maxsize=None)
def cap(x,y):
    h = gcd(x,y)
    a,b = x//h,y//h
    result = 1
    for p,q in factors(a*b):
        k,power = 0,1
        while power<q:
            k += 1
            power *= p
        result *= p**transport_cap(a,b,p,k)
    return result


def groups_for(d,c,height):
    forced = [q for p,q in factors(d) if c%q]
    for assignment in product(range(3),repeat=len(forced)):
        groups = tuple(prod(q for q,pos in zip(forced,assignment) if pos==i) for i in range(3))
        if groups[0]<=groups[1]<=groups[2]<=height:
            yield groups


def candidates(height):
    found = set()
    big = lcm(*range(1,height+1))
    for a,b in combinations_with_replacement(range(1,height+1),2):
        if gcd(a,b)>1:
            continue
        c = cap(a,b)
        for d in divisors(factors(lcm(c,big))):
            for groups in groups_for(d,c,height):
                modulus = prod(groups)
                steps = [modulus//v for v in groups]
                e = [k*pow(k,-1,v) if v>1 else 0 for k,v in zip(steps,groups)]
                for u in range(max(3,sum(steps)),(d-1)//(a+b)+1):
                    if gcd(u,d)>1:
                        continue
                    low = [v*u%modulus or modulus for v in e]
                    for x in range(low[0],u-low[1]-low[2]+1,modulus):
                        for y in range(low[1],u-x-low[2]+1,modulus):
                            z = u-x-y
                            if (z-low[2])%modulus:
                                continue
                            row = sorted((x,y,z))
                            if all(lcm(c,cap(i,j))>(a+b)*((i+j)//gcd(i,j))
                                   for i,j in combinations(row,2)):
                                found.add((a,b,d,u,*row))
    return found


@lru_cache(maxsize=8192)
def accepts(p,slopes,offsets,modulus,residue):
    """Read the fixed suffix, then all higher digits of a free multiplier."""
    def next_state(state,i):
        terms = tuple(t+s*i for s,t in zip(slopes,state))
        dest = tuple(x//p for x in terms)
        return dest if sum(terms)//p==sum(dest) else None

    state = offsets
    while modulus>1:
        residue,i = divmod(residue,p)
        state = next_state(state,i)
        if state is None:
            return False
        modulus //= p
    queue = [state]
    seen = {state}
    for state in queue:
        if state==(0,0,0):
            return True
        for i in range(p):
            dest = next_state(state,i)
            if dest is not None and dest not in seen:
                seen.add(dest)
                queue.append(dest)
    return False


def prime_possible(a,b,d,u,j,row,p):
    first = (a*u*j,b*u*j,d-(a+b)*u*j)
    zero = (0,0,0)
    one = (0,0,1)
    if u%p==0:
        w = dict(factors(u))[p]
        r = -pow(d,-1,w)%w
        return any(accepts(p,first,one,w*p,i) for i in range(r,w*p,w) if (d*i+1)%(w*p))
    if d%p==0:
        w = dict(factors(d))[p]
        for v,other in ((a,b),(b,a)):
            if other%w==0 and accepts(p,first,one,w,pow(v*u*j,-1,w)):
                return True
        if sum(x%w==0 for x in row)>=2:
            return accepts(p,row,zero,w,pow(u,-1,w))
        return False
    # For p | m, dividing the parts by its p-power after deleting the
    # isolated units digit gives precisely the homogeneous slope row.
    return any(accepts(p,s,o,p,t) for s,o in ((first,one),(first,zero),(row,zero))
               for t in range(1,p))


def run(height):
    shapes = candidates(height)
    print("height",height,"CRT shapes",len(shapes),flush=True)
    eliminated = Counter()
    survivors = []
    for a,b,d,u,*coords in sorted(shapes):
        for j in range(1,(d-1)//(u*(a+b))+1):
            if gcd(j,d)>1:
                continue
            for p in sorted({2,*dict(factors(d*u))}):
                if not prime_possible(a,b,d,u,j,tuple(coords),p):
                    eliminated[p] += 1
                    break
            else:
                survivors.append((a,b,d,u,j,*coords))
    print("excluded by prime",dict(sorted(eliminated.items())),"survivors",len(survivors),
          survivors[:12],flush=True)
    return shapes,eliminated,survivors


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--height",type=int,default=12)
    args = parser.parse_args()
    shapes,counts,remaining = run(args.height)
    if args.height==12:
        assert len(shapes)==27387
        assert remaining==[(4,11,13860,673,1,63,70,540),(7,9,2772,167,1,11,44,112)]
        for family in remaining:
            assert not any(r*x+s*y==t*z for x,y,z in permutations(family[-3:])
                           for r in range(1,13) for s in range(r,13) for t in range(1,13))
