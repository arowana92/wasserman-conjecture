"""Certify N=p**a*(p**j+d)+1 for j>=a>=1 and bounded positive d."""

from argparse import ArgumentParser
from collections import Counter
from functools import lru_cache
from math import gcd, isqrt

from check_close_powers import INITIAL_LIMIT, bases, check_complement, higher_weight, prime_base, subproducts
from check_six_bit_rays import divisors
from prime_certificates import primes_through
from search_even_endpoint import avoids, prime_powers


@lru_cache(maxsize=65536)
def factors_of(n):
    return tuple(divisors(n))


def sum_roots(limit):
    # N=Q*(Q-d)+1 vanishes mod U precisely when d=Q+Q**(-1) mod U.
    data = [{} for _ in range(limit+1)]
    for u in range(3, limit+1):
        for q in range(1, u):
            if gcd(q,u) == 1:
                data[u].setdefault((q+pow(q,-1,u)) % u, []).append(q)
    return data


def first_rows(P, Q, d, y, u, exceptional):
    # y is an unshifted small part after division by P. The unit is either on
    # the large part or the other small part. Both orientations are allowed.
    step = u*exceptional
    for x in range(step, d-y+1, step):
        yield (P*y, P*x, P*(Q-x-y)+1)
    # With the unit on x, x+y is divisible by exceptional and is Q mod u.
    start = exceptional*(Q*pow(exceptional,-1,u) % u)
    for t in range(y+1+(start-y-1) % step, d+1, step):
        yield (P*y, P*(t-y)+1, P*(Q-t))


def examine(P, p, Q, d, y, u, B, exceptional, counts, collector=None):
    n = P*Q+1
    T = Q//B
    if gcd(u,n//u) != 1 or gcd(T,B*P) != 1:
        return
    if any(2*w+1 > u for w in prime_powers_cached(T)):
        return
    counts["assignments"] += 1
    if collector is not None:
        collector.add((P,d,y,u,B,exceptional))
    required = (p,)+bases(u)+bases(B)
    for row in first_rows(P,Q,d,y,u,exceptional):
        counts["first rows"] += 1
        if all(avoids(n,row,r) for r in required):
            counts["first row passes"] += 1
            check_complement(n, u, T, counts)


@lru_cache(maxsize=32768)
def prime_powers_cached(n):
    return tuple(prime_powers(n).values())


def square(limit, minimum=INITIAL_LIMIT, collector=None):
    counts = Counter()
    roots = sum_roots(limit)
    for d in range(4, limit+1):
        for y in range(3, d):
            for u in factors_of(y):
                if u < 3:
                    continue
                residues = roots[u].get(d % u, ())
                max_t = isqrt(u**3//27)
                if not residues or (y//u)*d*max_t <= isqrt(minimum)-1:
                    continue
                # P is -d modulo Q. The factors exceptional at y divide d*y+1.
                for e in factors_of(d*y+1):
                    if e > d:
                        continue
                    for f in factors_of(y//u):
                        B = e*f
                        if gcd(B,u) > 1 or B*max_t <= isqrt(minimum)-1:
                            continue
                        inv = pow(B,-1,u)
                        for q in residues:
                            for T in range(q*inv % u,max_t+1,u):
                                Q = B*T
                                P = Q-d
                                n = P*Q+1
                                if P <= d or n <= minimum:
                                    continue
                                p = prime_base(P)
                                if p:
                                    examine(P,p,Q,d,y,u,B,e,counts,collector)
    return counts


def one_weight(P,p,H,d,counts,collector=None):
    # Start from N to check the higher-exponent branch. Enumerate the full
    # powers of N available to U and all unitary T dividing Q.
    Q = H+d
    n = P*Q+1
    counts["weights"] += 1
    for u in subproducts(w for w in prime_powers(n).values() if w < d):
        if not 3 <= u < d or 27*Q*Q > u*d**4:
            continue
        for T in subproducts(w for r,w in prime_powers(Q).items() if r != p and 2*w+1 <= u):
            B = Q//T
            if 27*T*T > u**3 or B*u >= d*d:
                continue
            for y in range(u,d,u):
                # Powers of B exceptional at this part divide P*y-1. All other
                # powers divide y/u. These two integers are coprime.
                e = gcd(B,P*y-1)
                if e <= d and (y//u) % (B//e) == 0:
                    examine(P,p,Q,d,y,u,B,e,counts,collector)


def higher(limit,minimum=INITIAL_LIMIT):
    counts = Counter()
    maximum = isqrt(limit**5//27)
    for p in primes_through(isqrt(maximum)):
        H = p*p
        while H < maximum:
            P = p
            while P < H:
                for d in range(4, min(limit,H-1)+1):
                    Q = H+d
                    if P*Q+1 > minimum and 27*Q*Q < d**5:
                        one_weight(P,p,H,d,counts)
                P *= p
            H *= p
    return counts


def preceding_multiples(limit):
    # Supplement the preceding-weight theorem where p divides the offset.
    counts = Counter()
    bound = limit**4//432
    for p in primes_through(limit):
        P = p
        while P <= bound:
            H = P
            while H*H <= bound*P:
                for d in range(p,min(limit,H-1)+1,p):
                    q = p
                    while d % (q*p) == 0:
                        q *= p
                    n = P*(H+d)
                    if n > INITIAL_LIMIT and 432*(n-1)**2*q <= P**3*d**4:
                        higher_weight(P*q,p,H//q,d//q,counts)
                H *= p
            P *= p
    return counts


if __name__ == "__main__":
    cli = ArgumentParser(description=__doc__)
    cli.add_argument("max_offset", nargs="?", type=int, default=1000)
    args = cli.parse_args()
    dmax = args.max_offset
    if dmax < 1 or 2*dmax*dmax+1 > INITIAL_LIMIT:
        cli.error("use a positive offset whose d>=H cases fit the initial theorem")
    print(f"Square case through {dmax}: {dict(square(dmax))}",flush=True)
    print(f"Higher powers through {dmax}: {dict(higher(dmax))}",flush=True)
    print(f"Preceding weights with p dividing d: {dict(preceding_multiples(dmax))}", flush=True)
