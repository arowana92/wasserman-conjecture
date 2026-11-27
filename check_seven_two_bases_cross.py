"""Cross-check the seven-factor theorem using integer log intervals and a joint modulus."""
from array import array
from bisect import bisect_left
from itertools import combinations
from math import gcd, lcm

from flint import arb, ctx

from check_seven_two_bases import simultaneous_approximants
from check_six_sevens_cross import ScaledLogs


def order(x, p):
    r, e = x % p, 1
    while r != 1:
        r = r*x % p
        e += 1
    return e


def residue_pairs(p, q):
    viable = []
    for a in range(lcm(6, order(p, q))):
        z = pow(p, a, q)
        if (6*z) % q+z >= q:
            continue
        for b in range(lcm(6, order(q, p))):
            w = pow(q, b, p)
            if pow(p, a, 504)*pow(q, b, 504) % 504 == 83 and (6*w) % p+w < p:
                viable.append((a, b))
    return viable


def masks(p, q):
    # Full sums are 6(1+t) with pq | t and t = 11 mod 72.
    m = 432*p*q
    n = 6*(1+p*q*(11*pow(p*q, -1, 72) % 72))
    seeds = [n % m]
    result = {}
    for length in range(6, 0, -1):
        present = bytearray(m)
        stack = array('I')
        for x in seeds:
            while not present[x]:
                stack.append(x)
                present[x] = 1
                x = 7*x % m
        result[length] = present
        seeds = ((z-1) % m for z in stack)
    return m, result


def certified_denominators(p, q, bound, logs):
    # LLL only proposes denominators. Verify each simultaneous approximation
    # against integer enclosures, then find a nonsingular submatrix directly.
    ctx.prec = 3*bound.bit_length()+510
    log7 = arb(7).log()
    candidates, _ = simultaneous_approximants(arb(p).log()/log7, arb(q).log()/log7, bound)
    intervals = [logs.base_seven(x) for x in (p, q)]
    good = []
    rows = []
    homogeneous = None
    scale = logs.scale
    for d in candidates:
        numerators = [(d*low+scale//2)//scale for low, _ in intervals]
        if any(128*(bound+3000)*max(abs(d*low-z*scale), abs(d*high-z*scale)) >= scale
               for z, (low, high) in zip(numerators, intervals)):
            continue
        good.append(d)
        x, y = numerators
        if homogeneous is None:
            for (u,v,w),(a,b,c) in combinations(rows, 2):
                if x*(v*c-w*b)-y*(u*c-w*a)+d*(u*b-v*a):
                    homogeneous = d
                    break
            rows.append((x,y,d))
    assert homogeneous is not None
    return good, homogeneous


def check_pass(p, q, module, masks_by_length, b, depth, next_b):
    bound = 10**b
    logs = ScaledLogs(3*bound.bit_length()+170)
    ds, homogeneous = certified_denominators(p, q, bound, logs)
    assert 16384*ds[-1] < 1 << 3000
    log6 = logs.base_seven(6)
    powers = [1]
    while powers[-1] < 16384*ds[-1]:
        powers.append(7*powers[-1])

    total = maximum = 0
    def search(j, r, gap):
        nonlocal total, maximum
        if not masks_by_length[j][r % module]:
            return
        total += 1
        reduced = r
        if r % 6 == 0:
            reduced //= 6
            for prime in (p, q, 7):
                while reduced % prime == 0:
                    reduced //= prime
        else:
            reduced = 0
        if reduced == 1:
            assert r.bit_length() < 3000
            d = homogeneous
        else:
            lo, hi = logs.base_seven(r)
            lo, hi = lo-log6[1], hi-log6[0]
            for d in ds:
                x, y = d*lo, d*hi
                k = (x//logs.scale)*logs.scale
                if min(x-k, k+logs.scale-y)*32 > logs.scale:
                    break
            else:
                raise ValueError(("increase precision", r))
        stop = bisect_left(powers, d*16384)
        if j == depth:
            maximum = max(stop, maximum)
            return
        for delta in range(gap, stop):
            search(j+1, r*powers[delta-gap]+1, delta)
    search(1, 1, 0)
    if next_b:
        assert (maximum+2)*(10**16*(2+3*b))**(6-depth) < 10**next_b
    return total


def small(p, q):
    # Generate exponent pairs independently, then compare the sum of binary
    # digit counts with the digit count of the weight (Legendre's formula).
    levels = [[base] for base in (p,q)]
    for sequence, base in zip(levels,(p,q)):
        while sequence[-1].bit_length() < 3000:
            sequence.append(sequence[-1]*base)
    count = 0
    for x in levels[0]:
        for y in levels[1]:
            t = x*y
            if (6*t).bit_length() > 3000:
                break
            if (t % 72, t % 7) == (11,6):
                assert bin(6*t).count('1')+bin(t+1).count('1') > bin(7*t+1).count('1')
                count += 1
    return count


def main():
    primes = [x for x in range(5, 48) if x != 7 and all(x%i for i in range(2, x))]
    positive = []
    for p,q in combinations(primes,2):
        if residue_pairs(p,q):
            positive.append((p,q))
    assert positive == [(17,43),(37,47),(41,43)]
    for p,q in positive:
        modulus, allowed = masks(p,q)
        counts = []
        for b,j,nb in ((130,2,82),(82,3,62),(62,4,42),(42,5,22),(22,6,0)):
            counts.append(check_pass(p,q,modulus,allowed,b,j,nb))
        print((p,q), "integer interval counts", counts, "finite pairs", small(p,q), flush=True)


if __name__ == "__main__":
    main()
