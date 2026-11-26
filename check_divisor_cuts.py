"""Reconstruct pairs of divisor-indexed binomials with small endpoint gcd."""
from argparse import ArgumentParser
from functools import lru_cache
from itertools import combinations_with_replacement
from math import comb, gcd, isqrt, lcm

from flint import fmpz

from search_even_endpoint import prime_powers


@lru_cache(maxsize=8192)
def divisors(n):
    return tuple(v for i in range(1, isqrt(n) + 1) if n % i == 0
                 for v in ((i,) if i*i == n else (i, n//i)))


@lru_cache(maxsize=32768)
def coefficient_gcd(n, a):
    value = 1
    for number in (n, n-1):
        for p, power in prime_powers(number).items():
            scale = p
            retained = 1
            while scale <= n and retained < power:
                if a % scale > n % scale:
                    retained *= p
                scale *= p
            value *= retained
    return value


def reconstruct(limit):
    hits = set()
    candidates = 0
    for common in range(2, limit+1):
        for k in range(1, common):
            for m in divisors(common-k):
                h = (common-k)//m
                for g in divisors(common):
                    kernel = h*(m*common-k*g)
                    if kernel <= 0:
                        continue
                    for x in divisors(kernel):
                        y = kernel//x
                        if x > y:
                            continue
                        s, rem = divmod(common+x, k)
                        u, other = divmod(common+y, k)
                        if rem or other or s < 2 or gcd(s,u) != g:
                            continue
                        n = m*(s//g)*u
                        b, a = n//u, n//s
                        if b < 2:
                            continue
                        assert n-1 < common**4
                        candidates += 1
                        if gcd(coefficient_gcd(n,a), coefficient_gcd(n,b)) == common:
                            hits.add((n,b,a,common))
    return candidates, hits


def direct(limit, bound):
    found = set()
    for n in range(4, bound+1):
        cuts = [(a, gcd(n*(n-1), comb(n,a))) for a in sorted(divisors(n)) if 1 < a < n]
        for (a,x), (b,y) in combinations_with_replacement(cuts,2):
            common = gcd(x,y)
            if common <= limit:
                found.add((n,a,b,common))
    return found


def rational_cuts(bound):
    tested = 0
    for n in range(2,bound+1):
        entries = []
        for a in range(1,n):
            m = gcd(n,a)
            r,s = a//m,n//m
            entries.append((r,s,r*(s-r),gcd(n*(n-1),comb(n,a))))
        for (r,s,f,x),(t,u,h,y) in combinations_with_replacement(entries,2):
            common = gcd(x,y)
            divisor = gcd(s,u)*((n-1)//gcd(n-1,lcm(f,h)))
            assert common % divisor == 0
            assert common*r*t > (n//lcm(s,u))*gcd(f,h)
            a,b = n*r//s,n*t//u
            if a < b:
                rr,tt = t,r
            else:
                rr,tt = r,t
            w = common*gcd(n-1,lcm(r,t))
            assert w > n//lcm(s,u)
            if min(a,b) >= 2:
                assert n-1 < 4*w**4*rr**2*tt
            tested += 1
    return tested


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('limit',type=int,nargs='?',default=128)
    parser.add_argument('--direct-through',type=int,default=2400)
    args=parser.parse_args()
    candidates, hits = reconstruct(args.limit)
    assert {h for h in hits if h[0] <= args.direct_through} == direct(args.limit,args.direct_through)
    print(f'Checked {rational_cuts(200):,} arbitrary cut pairs by exact binomial coefficients.')
    print(f'{candidates:,} reconstructed candidates, {len(hits):,} actual pairs for gcd <= {args.limit}.')
    if hits:
        n,a,b,common = max(hits)
        ambient = n*(n-1)
        assert fmpz.bin_uiui(n,a).gcd(ambient).gcd(fmpz.bin_uiui(n,b)) == common
    for bound in sorted({8,16,32,64,128,256,512,args.limit}):
        if bound > args.limit:
            continue
        limited = [h for h in hits if h[3] <= bound]
        if limited:
            print(bound,len(limited),max(limited))
