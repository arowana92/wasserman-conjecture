"""Exclude endpoint counterexamples with a bounded nonzero linear form on one row."""

import argparse
from collections import Counter
from functools import lru_cache
from itertools import combinations_with_replacement
from math import comb, gcd, isqrt, lcm

from search_endpoint import avoids, find_row


@lru_cache(maxsize=20000)
def factors(n):
    result = []
    p = 2
    while p*p<=n:
        if n%p==0:
            q=1
            while n%p==0:
                n//=p
                q*=p
            result.append((p,q))
        p = 3 if p==2 else p+2
    if n>1:
        result.append((n,n))
    return tuple(result)


def divisors(powers):
    result = [1]
    for p,full in powers:
        olds = result.copy()
        q = p
        while q<=full:
            result.extend(x*q for x in olds)
            q *= p
    return result


def first_row(normal,e,u,s_powers,s,v):
    r,b,neg_t = normal
    residues=[0]
    for p,q in s_powers:
        first = (e-r)%q==0
        other = (e-b)%q==0 or (e-neg_t)%q==0
        if not first and not other:
            return None
        if first:
            add = v*(s//q)*pow(s//q,-1,q)
            extra = [z+add for z in residues]
            residues = residues+extra if other else extra
    bases = [p for p,q in factors(u)+s_powers]
    for residue in residues:
        for x in range(residue%s or s,v-1,s):
            y, rem = divmod(e//u-neg_t*v-(r-neg_t)*x,b-neg_t)
            z = v-x-y
            if rem or min(y,z)<=0:
                continue
            row = (u*x,u*y,u*z)
            if all(avoids(row,p) for p in bases):
                return row
    return None


def scan(normal,radius,floor,row_finder=first_row):
    counts = Counter()
    for e in range(-radius,radius+1):
        if not e or abs(e)<=2:
            continue
        if e in normal:
            raise ValueError("radius includes a nonzero exceptional offset of magnitude above two")
        cap = lcm(*(abs(e-a) for a in normal))
        s_values = divisors(factors(cap))
        for u in divisors(factors(abs(e))):
            if u<3:
                continue
            limit = isqrt(u**3//27)
            for s in s_values:
                if s*limit+1<=floor or gcd(u,s)!=1:
                    continue
                start = max(1,(floor-1)//s+1)
                start += (-pow(s,-1,u)-start)%u
                for t in range(start,limit+1,u):
                    n = 1+s*t
                    v = n//u
                    if v<3 or gcd(u,v)>1 or gcd(s,t)>1:
                        continue
                    counts['assignments'] += 1
                    counts['max_weight'] = max(counts['max_weight'],n)
                    if any(2*q+1>u for p,q in factors(t)):
                        continue
                    counts['predecessor_power_test'] += 1
                    row = row_finder(normal,e,u,factors(s),s,v)
                    if row is None:
                        continue
                    counts['first_row_pass'] += 1
                    other = find_row(n,v,list(factors(v)),t,list(factors(t)))
                    if other is not None:
                        raise ArithmeticError((n,normal,e,row,other))
    return counts


def direct():
    checks = 0
    for n in range(3,30):
        rows = [(x,y,n-x-y,comb(n,x)*comb(n-x,y))
                for x in range(1,n-1) for y in range(1,n-x)]
        coefficients = {v for a,b,c,v in rows}
        for normal in [(1,2,-2),(2,3,-5),(-2,0,1),(1,1,1)]:
            for x,y,z,k in rows:
                e = sum(i*j for i,j in zip(normal,(x,y,z)))
                if not e or e in normal:
                    continue
                cap = lcm(*(abs(e-a) for a in normal))
                for h in coefficients:
                    g = gcd(n*(n-1),gcd(k,h))
                    assert g**3*abs(e)**3*cap**2>=27*(n-1)**2
                    checks += 1
    return checks


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--radius",type=int,default=2048)
    parser.add_argument("--floor",type=int,default=10**9,
                        help="omit previously certified weights up to this value")
    parser.add_argument("--small-check",action="store_true")
    args = parser.parse_args()
    if args.small_check:
        print("exact form/pair inequalities",direct())
    total_assignments = 0
    max_weight = 0
    for a,b,c in combinations_with_replacement(range(-2,3),3):
        if a==c or gcd(a,b,c)!=1 or (a,b,c)>(-c,-b,-a):
            continue
        normal = (b,c,a)
        counts = scan(normal,args.radius,args.floor)
        print(normal,dict(counts),flush=True)
        total_assignments += counts['assignments']
        max_weight = max(max_weight,counts['max_weight'])
    print(f"assignments={total_assignments}, largest weight={max_weight}")
