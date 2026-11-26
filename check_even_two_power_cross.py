"""Reconstruct the bounded cofactor assignments and test modular positives independently."""
from math import gcd, isqrt
from argparse import ArgumentParser
from itertools import groupby

from check_even_two_power_cofactors import positions, primes_and_powers, direct_progression


def least_factor(n):
    return next((p for p in range(2,isqrt(n)+1) if n%p == 0), n)


def reconstruct(D):
    w = D & -D
    if w < 8:
        return set()
    P = D//w
    if P == 1:
        return set()
    p = least_factor(P)
    k = P
    while k%p == 0:
        k //= p
    if k != 1:
        return set()
    shapes = set()
    for b in range(1,P):
        if gcd(P,b) > 1:
            continue
        for a in range(1,(D-w*b-1)//P+1,2):
            g = gcd(a,b)
            while g > 1:
                u = least_factor(g)
                while g%u == 0:
                    g //= u
                if P%u == 0 or a == w-1 == u:
                    continue
                shapes.add((w,P,u,a,b))
    return shapes


def audit(bound, small_multiplier_bound):
    largest = max(bound,small_multiplier_bound)
    primes,powers = primes_and_powers(largest//8)
    seen = set()
    count = 0
    w = 8
    while 5*w <= largest:
        limit = largest if w in (8,16) else bound
        for P, rows in groupby(positions(w,powers,primes,limit), key=lambda row: row[0]):
            generated = set()
            for _,p,u,a,b,r,modulus,slopes in rows:
                assert modulus == w*P*u and 0 < r < modulus
                assert a*P*r%w == b*w*r%P == 1 and (w*P*r+1)%u == 0
                assert min(slopes)>0 and sum(slopes)==w*P
                generated.add((w,P,u,a,b))
            assert generated == reconstruct(w*P)
            seen.add(w*P)
            count += len(generated)
        w *= 2
    for D in range(8,largest+1,8):
        if D not in seen and (D <= bound or (D & -D) in (8,16)):
            assert not reconstruct(D)
    print(count,"assignments independently reconstructed")

    # Exact factorial floors test the direct progression evaluator on finite
    # accepted examples including a large witness of the binary language.
    for base,slopes,r,step,k in ((2,(2,4,8),1,8,1),(3,(1,3,9),1,27,1),
                                 (2,(8043,896,3317),3,32,5)):
        q = r+step*k
        entries = (slopes[0]*q,slopes[1]*q,slopes[2]*q+1)
        layer = base
        while layer <= sum(entries):
            assert sum(entries)//layer == sum(x//layer for x in entries)
            layer *= base
        assert direct_progression(base,slopes,r,step)


if __name__ == '__main__':
    cli = ArgumentParser(description=__doc__)
    cli.add_argument('--bound',type=int,default=15423)
    cli.add_argument('--small-multiplier-bound',type=int,default=40000)
    args = cli.parse_args()
    audit(args.bound,args.small_multiplier_bound)
