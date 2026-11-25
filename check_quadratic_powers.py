"""Check the quadratic criterion at bounded coefficient pairs."""

from argparse import ArgumentParser
from collections import Counter
from itertools import combinations
from math import gcd, isqrt, lcm

from check_close_powers import INITIAL_LIMIT, bases, check_complement, prime_base, subproducts
from check_six_bit_rays import divisors
from search_even_endpoint import avoids, prime_powers


def shapes(h,d,eps):
    for a in range(h+1):
        for b in range(d+1):
            if a+b == 0:
                continue
            for c in range(h-a+1):
                for e in range(d-b+1):
                    F = ((a,b),(c,e),(h-a-c,d-b-e))
                    if F[1] < F[0] or not sum(F[2]) or (not eps and F[2] < F[1]):
                        continue
                    yield F


def factor_bounds(h,d,F,eps):
    D = [h*b-d*a for a,b in F]
    E = [abs(h*b*b-d*a*b+(2*eps-1)*a*a) if a else b for a,b in F]
    minors = {(i,j): F[i][0]*F[j][1]-F[i][1]*F[j][0] for i,j in combinations(range(3),2)}
    if eps == 0:
        v = gcd(*D,*minors.values())
        s = lcm(*(gcd(E[i],E[j],m) for (i,j),m in minors.items()))
        return v,s
    u = gcd(E[0],E[1],minors[0,1])
    groups = []
    for k in range(3):
        constants = [h*h*((i == 2)-(i == k))-d*D[i] for i in range(3)]
        constants += [D[i] for i in range(3) if (i != 2 and i != k) or k == 2]
        groups.append(gcd(*constants))
    return u,lcm(*groups)


def unshifted(h,d,F,vmax,smax,lo,counts):
    for v in divisors(vmax):
        residues_v = [x for x in range(v) if all((a*x+b) % v == 0 for a,b in F)]
        for s in divisors(smax):
            hi = (v**3*s*s-1)//(27*h*h)
            if hi < lo or gcd(v,s) != 1 or not residues_v:
                continue
            for y in range(s):
                if (y*(h*y+d)-1) % s:
                    continue
                parts = [a*y+b for a,b in F]
                allowed = lcm(*(gcd(parts[i],parts[j]) for i,j in combinations(range(3),2)))
                if allowed % s:
                    continue
                for x in residues_v:
                    start = y+s*((x-y)*pow(s,-1,v) % v)
                    for P in range(lo+(start-lo) % (s*v),hi+1,s*v):
                        counts["unshifted progression points"] += 1
                        p = prime_base(P)
                        if not p:
                            continue
                        Q = h*P+d
                        n = P*Q
                        if gcd(P*v,Q//v) != 1 or gcd(s,(n-1)//s) != 1:
                            continue
                        row = tuple(P*(a*P+b) for a,b in F)
                        if all(avoids(n,row,r) for r in (p,)+bases(v)+bases(s)):
                            counts["unshifted first rows"] += 1
                            check_complement(n,P*v,(n-1)//s,counts)


def shifted(h,d,F,umax,bmax,lo,hi,counts):
    for P in range(lo,hi+1):
        p = prime_base(P)
        if not p:
            continue
        n = P*(h*P+d)+1
        row = tuple(P*(a*P+b)+(i == 2) for i,(a,b) in enumerate(F))
        if not avoids(n,row,p):
            continue
        Q = h*P+d
        for u in divisors(gcd(umax,*row)):
            if u < 3 or gcd(u,n//u) != 1 or not all(avoids(n,row,r) for r in bases(u)):
                continue
            for T in subproducts(w for r,w in prime_powers(Q).items() if r != p and 2*w+1 <= u):
                B = Q//T
                if bmax % B or 27*T*T > u**3 or not all(avoids(n,row,r) for r in bases(B)):
                    continue
                counts["shifted first rows"] += 1
                check_complement(n,u,T,counts)


def check(limit):
    counts = Counter()
    for h in range(1,limit+1):
        for d in range(1,limit+1):
            for eps in (0,1):
                disc = d*d+(1-2*eps)*4*h
                square = disc >= 0 and isqrt(disc)**2 == disc
                universal = min(h,d) == 1 or (gcd(h,d) <= 2 and not square)
                special = {r for v in (h,d) if v > 1 for r in [prime_base(v)] if r}
                if not universal and not special:
                    continue
                if d*(h*d+d)+eps > INITIAL_LIMIT:
                    raise ValueError("initial interval must include P <= d")
                counts["all bases" if universal else "selected bases",eps] += 1
                lo = 1+(isqrt(d*d+4*h*(INITIAL_LIMIT-eps))-d)//(2*h)
                for F in shapes(h,d,eps):
                    if not universal:
                        upper,lower = zip(*F)
                        if not any(avoids(h,upper,r) and avoids(d,lower,r) for r in special):
                            continue
                    counts["shapes",eps] += 1
                    first,second = factor_bounds(h,d,F,eps)
                    if eps and 0 < first < 3:
                        continue
                    if min(first,second) < 1:
                        raise ArithmeticError(("vanishing bounds",h,d,F,eps))
                    if eps == 0:
                        hi = (first**3*second**2-1)//(27*h*h)
                    else:
                        hi = (isqrt(first**3*second**2//27)-d)//h
                    if hi < lo:
                        continue
                    counts["above initial interval",eps] += 1
                    if eps == 0:
                        unshifted(h,d,F,first,second,lo,counts)
                    else:
                        shifted(h,d,F,first,second,lo,hi,counts)
    print(dict(counts))


if __name__ == "__main__":
    cli = ArgumentParser(description=__doc__)
    cli.add_argument("max_coefficient",nargs="?",type=int,default=15)
    args = cli.parse_args()
    if args.max_coefficient < 1:
        cli.error("use a positive coefficient bound")
    check(args.max_coefficient)
