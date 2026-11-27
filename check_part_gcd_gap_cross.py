"""Independent scale-first enumeration of selected gaps a - gcd(b,c)."""

from argparse import ArgumentParser
from collections import Counter
from math import gcd, isqrt


def divisors(n):
    return [d for i in range(1,isqrt(n)+1) if n%i==0
            for d in ((i,n//i) if i*i!=n else (i,))]


def prime_powers(n):
    factors = []
    for p in range(2,isqrt(n)+1):
        if n%p==0:
            q = 1
            while n%p==0:
                q *= p
                n //= p
            factors.append((p,q))
    if n>1:
        factors.append((n,n))
    return factors


def cut_valuation(n,a,p):
    q = p
    v = 0
    while q<=n:
        v += n//q-a//q-(n-a)//q
        q *= p
    return v


def verify(h):
    if h in (0,1):
        raise ValueError("requires a nonzero gap different from one")
    counts = Counter()
    weights = set()
    H,E = abs(h),abs(h-1)
    for U in divisors(H):
        if U<3:
            continue
        tmax = isqrt(U**3//27)
        for g in divisors(E):
            max_m = (g*tmax*H + E//3)//(U*g)
            for m in range(1,max_m+1):
                M = U*g*m
                a = M+h
                if a<=0 or gcd(M,E)!=g:
                    continue
                # S divides g*a/U, and T S = N-1 >= 10^9.
                cap = g*(a//U)
                for k in range(1,cap*tmax//10**9 + 1):
                    S, rem = divmod(cap,k)
                    if rem:
                        continue
                    common = gcd(S,M)
                    if (h-1)%common:
                        continue
                    mod = M//common
                    residue = ((h-1)//common)*pow(S//common,-1,mod)%mod
                    low = (10**9+S-1)//S
                    for T in range(low+(residue-low)%mod,tmax+1,mod):
                        n = S*T+1
                        scale,rem = divmod(n-h,M)
                        if rem or scale<3 or gcd(U,n//U)>1 or gcd(S,T)>1:
                            continue
                        if any(2*q>=U for _,q in prime_powers(T)):
                            continue
                        if any(cut_valuation(n,a,p) for p,_ in prime_powers(U)):
                            continue
                        counts["possible assignments"] += 1
                        weights.add(n)
    return counts,weights


def exact_bounds(limit=70):
    """Check both quantitative bounds using exact trinomial coefficients and the minimum actual G."""
    from math import comb
    counts = Counter()
    for n in range(3,limit+1):
        rows = {(a,b,n-a-b): gcd(n*(n-1),comb(n,a)*comb(n-a,b))
                for a in range(1,n//3+1) for b in range(a,(n-a)//2+1)}
        divisors_seen = set(rows.values())
        minima = {v: min(gcd(v,w) for w in divisors_seen) for v in divisors_seen}
        for row, value in rows.items():
            G = minima[value]
            for i,a in enumerate(row):
                b,c = row[(i+1)%3],row[(i+2)%3]
                M = gcd(b,c)
                r = (b+c)//M
                for d in range(1,6):
                    for e in range(1,6):
                        h = d*a-e*M
                        if h in (0,d) or abs(h)>20:
                            continue
                        H,E = abs(h),abs(d-h)
                        if e+h*r == 0:
                            assert 432*H*(n-1)**2 <= G**3*E**2*e**4
                            counts["zero resultant"] += 1
                        else:
                            gap = 27*d*(n-1)-G**3*e*E**2*H**2*max(d,H)
                            assert gap<=0 or gap**2 <= 27*G**3*E**2*(E+H)**2*H
                            counts["nonzero resultant"] += 1
                            if gap>0:
                                counts["second term required"] += 1
    return counts


if __name__=="__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("gaps",type=int,nargs="+")
    args = parser.parse_args()
    print("Direct coefficient checks:",dict(exact_bounds()),flush=True)
    from check_binary_gap_finite import check_weight
    for h in args.gaps:
        counts,ws = verify(h)
        final = Counter()
        for n in sorted(ws):
            check_weight(n,final)
        print(h,dict(counts),"weights:",len(ws),"full endpoint checks:",dict(final),flush=True)
