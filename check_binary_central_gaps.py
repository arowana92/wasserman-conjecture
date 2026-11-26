"""Check the central equality for all power-of-two cut gaps in a binary weight range."""
from argparse import ArgumentParser
from collections import Counter
from math import comb, gcd, isqrt

from check_central_cuts import (binomial_carry, endpoint_primes, legendre, lower_powers_pass,
                                prime_list, quadratic_candidates, refinement_masks)


def candidates(bits, from_below=False, before_carries=False):
    # The two bounding quadratics are d^2 = t*N+s^2 and
    # (d+s)(d+s+1) = t*(N+d+s).
    s = 1
    while 16384*s**4 < (1 << bits)*(2*s+1):
        t = 2*s+1
        for m in range(max(2,(16384*s**4//t).bit_length()),bits+1):
            H = 1 << m
            if from_below:
                root = (isqrt(1+4*t*H)-1)//2
                first = root + (root*(root+1) < t*H) - s
                d = first + (3*s+1-first)%(8*s)
            else:
                last = isqrt(t*H+s*s-1)
                d = last - (last-3*s-1)%(8*s)
            if d <= s:
                continue
            n,rem = divmod(d*d-s*s,t)
            if rem or not H//2 <= n < H or H > n+d+s:
                continue
            a,rem = divmod(n-d-s,2)
            b = a+s
            if rem or not 0 < a < b <= n//2:
                continue
            assert 4*b*(n-a) == n*(n-1)
            if before_carries or ((a & (n-a)) == 0 and lower_powers_pass(n,a,b)):
                yield n,a,b
        s *= 2


def local_example():
    # Both rows satisfy every assigned lowest prime-power layer at this weight.
    # The higher binary carry is essential.
    s,i,j,m,f,g = 2,1,5,81,22,59
    v = i+8*s*m
    u = (j*v+2*s)//i
    Sodd, T = (v+i)//2,(u-j)//(2*s)
    S = 4*s*Sodd
    n = u*v
    A = (4*s*u*f,4*s*u*g,u*Sodd)
    B = (v*T,v*T,n-2*v*T)
    a,b = sum(A[:2]),sum(B[:2])
    assert S*T == n-1 and gcd(u,v) == gcd(S,T) == 1
    assert sum(A) == sum(B) == n and b-a == s and 4*b*(n-a) == n*(n-1)
    for p in prime_list(S):
        for factor,row in ((S,A),(T,B)):
            if factor%p:
                continue
            power = p
            while factor%(p*power) == 0:
                power *= p
            assert sorted(x%power for x in row) == [0,0,1]
    assert all(x%u == 0 for x in A) and all(x%v == 0 for x in B)
    assert legendre(n,2) > sum(legendre(x,2) for x in A)
    return n


def lifting_controls():
    d,s = 2,1
    for _ in range(100):
        d,s = 649*d+2340*s,180*d+649*s
        n = 4*s*s-3
        assert d*d == 13*s*s-9
        for x in (d-s,d+s):
            a = (n-x)//2
            assert a & (n-a)


def odd_gap_candidates(bits,odd_limit):
    for odd in range(1,odd_limit+1,2):
        q = 1
        while 2*q*odd < 1 << bits:
            s = q*odd
            t = 2*s+1
            step = 8*q
            for m in range((2*s+1).bit_length(),bits+1):
                H = 1 << m
                first = isqrt(t*H)-s
                last = isqrt(t*H+s*s-1)
                for sign in (-1,1):
                    d0 = first+(sign*(3*s+1)-first)%step
                    for d in range(d0,last+1,step):
                        n,rem = divmod(d*d-s*s,t)
                        if rem or not H//2 <= n < H or H > n+d+s:
                            continue
                        a,rem = divmod(n-d-s,2)
                        b = a+s
                        if rem or not 0 < a < b <= n//2:
                            continue
                        if (a & (n-a)) and (b & (n-b)):
                            continue
                        if lower_powers_pass(n,a,b):
                            assert 4*b*(n-a) == n*(n-1)
                            yield n,a,b
            q *= 2


def odd_gap_controls(bits,odd_limit):
    cases = set(odd_gap_candidates(bits,odd_limit))
    check_bits = min(128,bits)
    control = {case for case in quadratic_candidates(512,check_bits)
               if (case[2]-case[1])//((case[2]-case[1]) & -(case[2]-case[1])) <= odd_limit
               and lower_powers_pass(*case)}
    assert control == {case for case in cases if case[0].bit_length() <= check_bits and
                       case[2]-case[1] <= 512}
    primes = prime_list(10000)
    exceptions = []
    for n,a,b in sorted(cases):
        for p in primes:
            if n%p <= 1 and binomial_carry(n,a,p) and binomial_carry(n,b,p):
                assert all(legendre(n,p) > legendre(z,p)+legendre(n-z,p) for z in (a,b))
                break
        else:
            if n > 1_000_000:
                raise ArithmeticError(f"Unresolved odd-part case at {n.bit_length()} bits.")
            common = gcd(n*(n-1),gcd(comb(n,a),comb(n,b)))
            if common > 1:
                continue
            ps = endpoint_primes(n)
            left,right = (refinement_masks(n,z,ps) for z in (a,b))
            assert all(u & v for u in left for v in right)
            exceptions.append((n,a,b))
    print(f"Odd gap parts through {odd_limit}, weights below 2**{bits}: "
          f"{len(cases)} candidates, binomial exceptions {exceptions}.")


def verify(bits,cross_bits):
    full = set(candidates(bits))
    cross = {case for case in full if case[0].bit_length() <= cross_bits}
    assert cross == set(candidates(cross_bits,from_below=True))
    geometric = set(candidates(cross_bits,before_carries=True))
    assert geometric == set(candidates(cross_bits,from_below=True,before_carries=True))
    witnesses = Counter()
    primes = prime_list(10000)
    for n,a,b in sorted(full):
        for p in primes:
            if n%p <= 1 and binomial_carry(n,a,p) and binomial_carry(n,b,p):
                assert all(legendre(n,p) > legendre(z,p)+legendre(n-z,p) for z in (a,b))
                witnesses[p] += 1
                break
        else:
            raise ArithmeticError(f"Unresolved equality at {n.bit_length()} bits and gap {b-a}.")
    print(f"All power-of-two equality gaps below weight 2**{bits}, {len(full)} candidates.")
    print(f"Second enumeration below 2**{cross_bits}: {len(geometric)} before carry "
          f"filters and {len(cross)} after them.")
    print(f"Common primes with factorial-valuation checks: {dict(witnesses)}.")


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--bits",type=int,default=4096)
    parser.add_argument("--cross-bits",type=int,default=512)
    parser.add_argument("--odd-bits",type=int,default=512)
    parser.add_argument("--odd-limit",type=int,default=15)
    args = parser.parse_args()
    if not 2 <= args.cross_bits <= args.bits or args.odd_bits < 2 or args.odd_limit < 1:
        parser.error("Require 2 <= cross-bits <= bits, odd-bits >= 2 and odd-limit >= 1.")
    print(f"Lowest-layer positive control at weight {local_example()}.")
    lifting_controls()
    verify(args.bits,args.cross_bits)
    odd_gap_controls(args.odd_bits,args.odd_limit)
