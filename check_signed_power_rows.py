"""Test positive and negative single-part shifts of scaled prime-power triples."""

from argparse import ArgumentParser
from collections import Counter
from functools import lru_cache
from math import gcd, isqrt


INITIAL = 10**9


@lru_cache(None)
def divisors(n):
    small = [i for i in range(1,isqrt(n)+1) if n%i == 0]
    return small + [n//i for i in reversed(small) if i*i != n]


@lru_cache(None)
def prime_divisors(n):
    primes = []
    for q in range(2,isqrt(n)+1):
        if n%q == 0:
            primes.append(q)
            while n%q == 0:
                n //= q
    if n > 1:
        primes.append(n)
    return primes


def powers_through(limit):
    sieve = bytearray(b"\x01")*(limit+1)
    for p in range(2,isqrt(limit)+1):
        if sieve[p]:
            sieve[p*p::p] = bytes((limit - p*p)//p + 1)
    result = []
    for p in range(2,limit+1):
        if sieve[p]:
            q = p
            while q <= limit:
                result.append((q,p))
                q *= p
    return sorted(result)


def avoids(row, p):
    a,b,c = row
    while a or b or c:
        if a%p + b%p + c%p >= p:
            return False
        a,b,c = a//p,b//p,c//p
    return True


def weight_prime(h, powers, counts):
    """The prime of the direction divides N. A large U forces the shift at its only unit."""
    H,B = abs(h),abs(h-1)
    if H**3*B**4 < 728*INITIAL:
        return set()
    assert (16*H)**3 * B**2 < 27*INITIAL**2
    pending = set()
    for q,p in powers:
        if 6075*q >= 64*B*B*H:
            break
        # Enumerate all divisors of H*q with nontrivial p-part.
        if (H*q)**3*B**2 < 27*INITIAL**2:
            continue
        qmax = isqrt(256*B*B*H*q//6075)+1
        hfree = H
        while hfree%p == 0:
            hfree //= p
        for u in divisors(hfree):
            pp = p
            while pp <= H*q//hfree:
                U = u*pp
                pp *= p
                if U**3*B*B < 27*INITIAL**2:
                    continue
                hi_n = 1+isqrt(U**3*B*B//27)
                r = q
                while r <= qmax:
                    s = 1+q+r
                    # n = s*U*k - (s-1)*h with positive shifted part U*k.
                    correction = (s-1)*h
                    lo_k = max(1,(INITIAL+correction)//(s*U)+1)
                    hi_k = (hi_n+correction)//(s*U)
                    for k in range(lo_k,hi_k+1):
                        M = U*k-h
                        n = s*M+h
                        g = gcd(M,B)
                        if U**3*g*g < 27*(n-1)**2 or gcd(U,n//U)>1:
                            continue
                        counts["weight-prime size candidates"] += 1
                        row = (M+h,M*q,M*r)
                        if all(avoids(row,d) for d in prime_divisors(u)) and avoids(row,p):
                            pending.add(tuple(sorted(row)))
                    r *= p
    return pending


def other_prime(h, counts):
    """When p does not divide N, enumerate the bounded normalized scale and the low residues."""
    H,B = abs(h),abs(h-1)
    assert B*B < INITIAL and B*isqrt(H**3//27) < INITIAL
    pending = set()
    for U in divisors(H):
        if U<3:
            continue
        tmax = isqrt(U**3//27)
        for m in range(U,tmax,U):
            for g in divisors(B):
                beta = (h-1)//g
                if gcd(m,beta)>1:
                    continue
                pcap = (m+1)*(m+abs(beta))
                if g*tmax*pcap < INITIAL:
                    continue
                primes = set()
                for k in (1,2,3):
                    primes.update(prime_divisors(abs(k*m+beta)))
                M = m*g
                rmax = (tmax*pcap + abs(beta))//m
                for p in primes:
                    r = 1
                    while r <= rmax:
                        q = 1
                        while q <= r:
                            s = 1+q+r
                            n = s*M+h
                            if n>INITIAL and n%p and gcd(n//U,U)==1:
                                cap = g*gcd(m*s+beta,r)
                                if U**3*cap**2 >= 27*(n-1)**2:
                                    counts["other-prime size candidates"] += 1
                                    for slot in range(3):
                                        row = [M*z+int(i==slot)*h for i,z in enumerate((1,q,r))]
                                        if min(row)>0 and all(avoids(row,d) for d in prime_divisors(U)):
                                            pending.add(tuple(sorted(row)))
                            q *= p
                        r *= p
    return pending


def run(limit):
    if not 2 <= limit <= 512:
        raise ValueError("use a shift bound in 2..512")
    powers = powers_through(64*limit*(limit+1)**2//6075)
    totals = Counter()
    possible = set()
    for h in range(-limit,limit+1):
        if h in (0,1):
            continue
        left = weight_prime(h,powers,totals)
        right = other_prime(h,totals)
        totals["weight-prime distinct weights per shift"] += len({sum(row) for row in left})
        totals["other-prime distinct weights per shift"] += len({sum(row) for row in right})
        possible.update(left|right)
    return totals,possible


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--limit",type=int,default=128)
    args = parser.parse_args()
    counts, rows = run(args.limit)
    weights = {sum(row) for row in rows}
    print(dict(counts), "distinct weights:",len(weights),"maximum:",max(weights,default=0),flush=True)
    if weights:
        from check_binary_gap_finite import check_weight
        final = Counter()
        for weight in sorted(weights):
            check_weight(weight,final)
        print("Endpoint tests completed:",dict(final))
        from check_signed_power_companions import check_rows
        print("Independent check on",len(rows),"rows:",dict(check_rows(rows)))
