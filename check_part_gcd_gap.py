"""Bound rational approximations between one part and the gcd of the other two."""

from argparse import ArgumentParser
from collections import Counter
from functools import lru_cache
from math import gcd, isqrt, lcm

from check_signed_power_rows import INITIAL, divisors, prime_divisors


@lru_cache(None)
def predecessor_products(u):
    """Divisors coprime to u supported on two entries of some positive partition."""
    caps = set()
    for a in range(1,u//3+1):
        for b in range(a,(u-a)//2+1):
            c = u-a-b
            cap = lcm(gcd(a,b),gcd(a,c),gcd(b,c))
            while (common := gcd(cap,u)) > 1:
                cap //= common
            caps.add(cap)
    return sorted({t for cap in caps for t in divisors(cap)})


def cut_avoids(part, n, p):
    while part:
        if part%p > n%p:
            return False
        part,n = part//p,n//p
    return True


def k_bound(g, t, u, d, e, h):
    H,E = abs(h),abs(d-h)
    max_scale = g*t*max(d,H) + E//(e+2*d)
    return g*t*(e*max_scale+H)//(d*u*INITIAL)


def scan(h, d, e, counts):
    H,E = abs(h),abs(d-h)
    weights = set()
    for u in divisors(H):
        if u<3:
            continue
        tmax = isqrt(u**3//27)
        # At delta = e + h*r = 0, r=x+y is fixed and S <= E*x*y.
        if h<0 and e%H==0 and e//H>=2:
            if E*tmax*(e//H)**2 >= 4*INITIAL:
                raise ValueError("zero resultant exceeds the initial range at these parameters")
        if not k_bound(E,tmax,u,d,e,h):
            continue
        for t in predecessor_products(u):
            for g in divisors(E):
                beta = (h-d)//g
                for k in range(1,k_bound(g,t,u,d,e,h)+1):
                    remainder = (h//u)*t-k*beta
                    if not remainder:
                        continue
                    for m in divisors(abs(remainder)):
                        slope, rem = divmod(e*g*t + remainder//m,k*u)
                        if rem:
                            continue
                        r, rem = divmod(slope-e,d)
                        if rem or r<2:
                            continue
                        M = g*u*m
                        a, rem = divmod(e*M+h,d)
                        n = a + M*r
                        if rem or a<1 or n<=INITIAL or gcd(M,a-1)!=g:
                            continue
                        assert (n-1)*k*u == g*t*a
                        if a%u or (n-1)%t or gcd(t,(n-1)//t)>1 or gcd(u,n//u)>1:
                            continue
                        counts["assignments after unitary tests"] += 1
                        if all(cut_avoids(a,n,p) for p in prime_divisors(u)):
                            counts["assignments after binomial cuts"] += 1
                            weights.add(n)
    return weights


def range_check(limit, denominators, numerators):
    counts = Counter()
    weights = set()
    for d in range(1,denominators+1):
        for e in range(1,numerators+1):
            if gcd(e,d)>1:
                continue
            for h in range(-limit,limit+1):
                if h not in (0,d):
                    weights.update(scan(h,d,e,counts))
    return counts, weights


def analytic_table():
    for bound, H in ((1,121),(8,80),(16,69),(32,60),(64,52),(100,48),
                     (128,45),(256,35),(512,13),(1000,5)):
        # Bound the two summands A+B*sqrt(H/27) by exact squaring, for each d.
        # e <= bound. Use H>=|h| and E <= d+H at every step.
        for d in range(1,bound+1):
            E = d+H
            delta = 27*d*INITIAL-bound*E**2*H**2*max(d,H)
            assert delta>0 and delta**2 > 27*H*E**2*(E+H)**2
            assert E*bound**2 < 20*INITIAL  # 4*sqrt(27)>20
        print("analytic bound",bound,"remainder through",H)


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--limit",type=int,default=512)
    parser.add_argument("--denominators",type=int,default=1)
    parser.add_argument("--numerators",type=int,default=1)
    args = parser.parse_args()
    analytic_table()
    counts, weights = range_check(args.limit,args.denominators,args.numerators)
    print(dict(counts),"weights:",len(weights),"maximum:",max(weights,default=0),flush=True)
    if weights:
        from check_binary_gap_finite import check_weight
        final = Counter()
        for n in sorted(weights):
            check_weight(n,final)
        print("Full endpoint tests:",dict(final))
