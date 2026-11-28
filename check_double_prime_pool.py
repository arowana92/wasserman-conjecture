"""Produce disjoint possible-common-prime sets at weights N=2p."""

import argparse
from array import array
from math import comb, gcd

from check_consecutive_certificate_lower_bound import primes_through
from check_direction_pair_count import valuation


def avoiding_row(n, p):
    # Separate two digit units, leaving at least one in the third part.
    rest = n
    first = []
    for _ in range(2):
        power = 1
        while rest % (p*power) == 0:
            power *= p
        first.append(power)
        rest -= power
        if rest == 0:
            return None
    return (*first,rest)


def finite(bound):
    primes = primes_through(bound)
    least = array("I",[0])*(bound+1)
    for p in primes:
        if p*p > bound:
            break
        for n in range(p*p,bound+1,p):
            if not least[n]:
                least[n] = p
    used = set()
    witnesses = seen = exact = 0
    for p in reversed(primes):
        n = 2*p
        if not bound//2 < n <= bound or p == 2:
            continue
        odd = avoiding_row(n,2)
        if odd is None:
            continue
        m = n-1
        r = least[m]
        q = m//r if r else m
        if least[q]:
            continue
        row = avoiding_row(n,r) if r and r != q else odd
        assert row is not None and valuation(n,odd,2) == 0
        if r and r != q:
            assert valuation(n,row,r) == 0
        seen += 1
        if q in used or p in used:
            continue
        used.update((p,q))
        witnesses += 1
        assert q*q >= m and p*p > m
        # The endpoint contains every candidate common prime. Evaluate their
        # powers directly, including a repeated prime factor of m.
        common = 1
        for t in {2,p,q} | ({r} if r else set()):
            endpoint = n*m
            power = 1
            while endpoint%t == 0:
                power *= t
                endpoint //= t
            common *= gcd(power,t**min(valuation(n,odd,t),valuation(n,row,t)))
        assert common%p == 0 and p*q*q % common == 0
        if n <= 4000:
            def coefficient(triple):
                a,b,_ = triple
                return comb(n,a)*comb(n-a,b)
            assert gcd(n*m,coefficient(odd),coefficient(row)) == common
            exact += 1
    assert len(used) == 2*witnesses
    print("weight bound",bound,"eligible weights in upper half",seen,
          "disjoint prime pairs",witnesses,"exact coefficient checks",exact)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bound",type=int,default=10000000)
    args = parser.parse_args()
    finite(4000)
    finite(args.bound)
