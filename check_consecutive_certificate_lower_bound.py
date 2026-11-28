"""Construct consecutive weights whose three explicit rows need many different primes."""

import argparse
from math import comb, gcd, isqrt, prod

from check_direction_pair_count import valuation


def primes_through(limit):
    sieve = bytearray(b'\x01')*(limit+1)
    sieve[:2] = b'\x00\x00'
    for p in range(2,isqrt(limit)+1):
        if sieve[p]:
            sieve[p*p::p] = bytes((limit-p*p)//p+1)
    return [i for i in range(2,limit+1) if sieve[i]]


def next_power(p, value):
    q = 1
    while q <= value:
        q *= p
    return q


def construct(length, cutoff, start):
    high = length+1
    if cutoff < high or length < 1:
        raise ValueError("require cutoff >= length+1 >= 2")
    primes = primes_through(cutoff)
    low = {p:next_power(p,high) for p in primes if p <= high}
    needed = dict(low)
    parts = {}
    for h in range(2,high+1):
        active = [p for p in low if h*(h-1) % p == 0]
        groups = ([p for p in active if p**3 <= high**2],
                  [p for p in active if p**3 > high**2])
        pair = []
        for group in groups:
            a = prod(low[p] for p in group)
            pair.append(a)
            for p in group:
                needed[p] = max(needed[p],a)
        parts[h] = pair
    residue, modulus = 0, 1
    for p in primes:
        q = next_power(p,needed[p]) if p in low else p
        target = -low[p] if p in low else 0
        residue += modulus*((target-residue)*pow(modulus,-1,q) % q)
        modulus *= q
    start = max(start,1+max(a for pair in parts.values() for a in pair))
    return start+(residue-start)%modulus,modulus,parts,primes


def verify(length, cutoff, start, above_power=False):
    base,modulus,parts,primes = construct(length,cutoff,start)
    primorial = prod(primes)
    exact = 0
    quotients = []
    for h,(a,b) in parts.items():
        n = base+h
        rows = ((1,1,n-2),(h,a,base-a),(h,b,base-b))
        assert all(min(row) > 0 and sum(row) == n for row in rows)
        assert n%h == (n-1)%(h-1) == 0
        q = n*(n-1)//(h*(h-1))
        assert gcd(q,primorial) == 1
        for p in primes:
            if n*(n-1)%p == 0:
                assert any(valuation(n,row,p) == 0 for row in rows[1:])
        if n <= 30000:
            assert gcd(*(comb(n,x)*comb(n-x,y) for x,y,z in rows)) == q
            exact += 1
        if above_power:
            assert 0 <= n-start < start and (n-start)**7 <= 6**6*(n-1)**4
        for prior_h,prior in enumerate(quotients,2):
            expected = (n-1)//(h-1) if prior_h+1 == h else 1
            assert gcd(prior,q) == expected
        quotients.append(q)
    return base,modulus,exact


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--length",type=int,default=128)
    parser.add_argument("--cutoff",type=int,default=1000)
    parser.add_argument("--exponent",type=int,default=4096)
    args = parser.parse_args()
    exact = 0
    for h in range(1,7):
        for y in range(h+1,17):
            for x in (1,11,101):
                _,_,count = verify(h,y,x)
                exact += count
    start = 2**args.exponent
    base,modulus,_ = verify(args.length,args.cutoff,start,True)
    print("exact small coefficient gcds",exact)
    print("length",args.length,"avoided prime cutoff",args.cutoff,"binary exponent",args.exponent)
    print("CRT modulus bits",modulus.bit_length(),"offset bits",(base-start).bit_length())
    print("minimum distinct primes for the explicit witnesses",(args.length+1)//2)
