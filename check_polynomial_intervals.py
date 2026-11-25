"""Close the finite interval checks for N=P^4+A*P^3+P^2+P+c at dyadic A and P."""

from argparse import ArgumentParser
from collections import Counter
from math import gcd, isqrt, lcm

from check_close_powers import INITIAL_LIMIT, check_complement
from search_even_endpoint import avoids, prime_powers


def submasks(n):
    m = n
    while m:
        yield m
        m = (m-1) & n


def check_weight(n, counts):
    counts["weights"] += 1
    m = n-1
    target = 27*m*m
    possible = set()
    for z in submasks(n):
        if not 3 <= z <= n-6:
            continue
        u = gcd(n,z)
        if u < 3:
            continue
        s = gcd(m,z*(z-1))
        if u**3*s*s >= target:
            possible.add(z)
    counts["eligible single parts"] += len(possible)
    for a in sorted(possible):
        if 3*a > n:
            break
        for b in submasks(n^a):
            c = n-a-b
            if a >= b or b >= c or b not in possible or c not in possible:
                continue
            counts["binary rows"] += 1
            u = gcd(a,b,c)
            s = gcd(m,lcm(gcd(a,b),gcd(b,c),gcd(c,a)))
            if u**3*s*s < target:
                continue
            row = (a,b,c)
            assigned_u = assigned_s = 1
            for base,w in prime_powers(n).items():
                if u % w == 0 and avoids(n,row,base):
                    assigned_u *= w
            for base,w in prime_powers(m).items():
                if s % w == 0 and avoids(n,row,base):
                    assigned_s *= w
            if assigned_u < 3 or n//assigned_u < 3:
                continue
            T = m//assigned_s
            if 27*T*T <= assigned_u**3:
                counts["complement tests"] += 1
                check_complement(n,assigned_u,T,counts)


def check(A):
    # Resultants for bulk offsets use coefficient l2 norm at most this value.
    norm = isqrt(A*A+3+(A-4)**2)+1
    counts = Counter()
    P = 2
    while P <= norm**6:
        center = P**4+A*P**3+P**2+P
        for c in range(2,A-3):
            if center+c > INITIAL_LIMIT:
                check_weight(center+c,counts)
        P *= 2
    print(f"A={A}, analytic cutoff {norm}^6, last tested exponent {(P//2).bit_length()-1}: {dict(counts)}")


if __name__ == "__main__":
    cli = ArgumentParser(description=__doc__)
    cli.add_argument("coefficient",type=int,nargs="?",default=1024)
    args = cli.parse_args()
    A = args.coefficient
    if A < 8 or A & (A-1):
        cli.error("use a power of two of at least 8")
    check(A)
