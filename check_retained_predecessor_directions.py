"""Test the full missing predecessor divisor, including partially lost powers."""

import argparse
from bisect import bisect_right
from functools import lru_cache
from itertools import combinations
from math import gcd,lcm,log

from carry_blocks import CarryGraph
from check_affine_endpoint import factors
from check_direction_pair_count import divisors,heavy_rows
from check_direction_predecessor_cross import raw_accepts
from check_direction_predecessor_loss import prime_loss,vp


@lru_cache(maxsize=50000)
def profile(row):
    h=lcm(*(gcd(a,b) for a,b in combinations(row,2)))
    return dict(factors(h))


def group_cases(limit):
    cases=[]
    total=0
    for row,t in heavy_rows(limit):
        total+=1
        H=profile(row)
        if all(prime_loss(row,p,H[p])[0]%w==0 for p,w in factors(t)):
            cases.append((*row,t))
    return total,sorted(cases)


def direct_cases(limit):
    out=[]
    for x in range(1,limit//3+1):
        for y in range(x,(limit-x)//2+1):
            d=gcd(x,y)
            if d==1:
                continue
            for z in range(y,limit-x-y+1):
                if gcd(d,z)>1:
                    continue
                u=x+y+z
                h=d*gcd(x,z)*gcd(y,z)
                if h<=u:
                    continue
                row=(x,y,z)
                E=h
                for p,w in factors(h):
                    e=vp(w,p)
                    for j in range(e,0,-1):
                        if raw_accepts(row,p,e,j):
                            break
                        E//=p
                    if E<=u:
                        break
                if E>u:
                    out.extend((*row,t) for t in divisors(E) if t>u)
    return sorted(out)


def prime_tests(row,t):
    H=profile(row)
    result=[]
    for p,q in factors(t):
        e=vp(H[p],p)
        unit,=[i for i,x in enumerate(row) if x%p]
        low=next(i for i,x in enumerate(row) if x%p==0 and vp(x,p)==e)
        off=tuple(int(i in (unit,low)) for i in range(3))
        result.append((p,q,e-vp(q,p),off))
    return result


def path_count(row,p,q,budget,off,length):
    r=pow(sum(row),-1,q)
    start=tuple((r*x-d)//q for x,d in zip(row,off))
    if budget==0:
        return CarryGraph(p,row,start).block_bound(length)
    states=[(start,0,budget)]
    indices={states[0]:0}
    edges=[]
    for co,incoming,left in states:
        next_states=[]
        for digit in range(p):
            vs=[v+a*digit for a,v in zip(row,co)]
            carry=(sum(v%p for v in vs)+incoming)//p
            if carry>left:
                continue
            dest=(tuple(v//p for v in vs),carry,left-carry)
            if dest not in indices:
                indices[dest]=len(states)
                states.append(dest)
            next_states.append(indices[dest])
        edges.append(next_states)
    paths=[1]*len(states)
    for _ in range(length):
        paths=[sum(paths[i] for i in dests) for dests in edges]
    return sum(paths)


def scan(row,t,bound,reverse):
    u=sum(row)
    r=pow(u,-1,t)
    tests=[]
    for p,q,budget,off in prime_tests(row,t):
        levels=[p]
        while levels[-1]<bound:
            levels.append(levels[-1]*p)
        tests.append((levels,budget,off))
    last=bound//u
    v=last-(last-r)%t if reverse else r
    count=0
    points=[]
    while 0<v<=last:
        n=u*v-2
        parts=[tuple(x*v-d for x,d in zip(row,off)) for _,_,off in tests]
        lengths=[bisect_right(seq,n) for seq,_,_ in tests]
        costs=[0]*len(tests)
        ends=[None]*len(tests)
        jump=None
        for depth in range(1,max(lengths,default=0)+1):
            for i,((seq,budget,off),k) in enumerate(zip(tests,lengths)):
                if depth>k:
                    continue
                q=seq[k-depth]
                top=n//q
                floors=[x//q for x in parts[i]]
                cost=top-sum(floors)
                if not cost:
                    continue
                costs[i]+=cost
                if reverse:
                    first=(top*q+2+u-1)//u
                    ends[i]=first if ends[i] is None else max(ends[i],first)
                    if costs[i]>budget:
                        jump=ends[i] if jump is None else min(jump,ends[i])
                else:
                    stop=min(((f+1)*q+d+x-1)//x for x,d,f in zip(row,off,floors))
                    ends[i]=stop if ends[i] is None else min(ends[i],stop)
                    if costs[i]>budget:
                        jump=ends[i] if jump is None else max(jump,ends[i])
            if jump is not None:
                break
        count+=1
        if jump is None:
            points.append(v)
            jump=v if reverse else v+1
        v=jump-1-(jump-1-r)%t if reverse else jump+(r-jump)%t
    return count,sorted(points)


def run(limit,digits,independent):
    total,cases=group_cases(limit)
    print("sum bound",limit,"coarse assignments",total,"retained families",len(cases),
          "directions",len({c[:3] for c in cases}),flush=True)
    if independent:
        assert cases==direct_cases(limit)
        print("coordinate and raw-valuation enumeration agrees",flush=True)
    maximum=forward=backward=0
    bound=10**digits
    for *parts,t in cases:
        row=tuple(parts)
        nums=[]
        for p,q,b,off in prime_tests(row,t)[:2]:
            paths=path_count(row,p,q,b,off,96)
            k=int(1000*log(paths,p)/96)+1
            assert paths**1000<p**(96*k)
            nums.append(k)
        maximum=max(maximum,sum(nums))
        assert len(nums)==2 and sum(nums)<1000
        l,one=scan(row,t,bound,False)
        r,two=scan(row,t,bound,True)
        assert one==two
        if one:
            raise ArithmeticError((row,t,one))
        forward+=l
        backward+=r
    print("worst dimension sum numerator",maximum,"decimal exponent",digits,
          "forward/backward intervals",forward,backward,flush=True)
    if limit==2000:
        assert total==133318 and len(cases)==2243 and maximum==966


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--limit",type=int,default=2000)
    parser.add_argument("--digits",type=int,default=1000)
    parser.add_argument("--independent",action="store_true")
    args=parser.parse_args()
    run(args.limit,args.digits,args.independent)
