"""Exact minima near N/t for binomial cuts against arbitrary positive trinomials."""

from collections import Counter
from functools import lru_cache
from math import comb, gcd, isqrt, lcm

from flint import fmpz

from carry_blocks import CarryGraph, intersection
from check_central_cuts import binomial_carry
from search_endpoint import find_row


RADII = {2: 4616, 3: 4048, 4: 113, 5: 2912, 6: 2344}


@lru_cache(maxsize=100000)
def factors(n):
    return tuple((int(p), int(p)**int(e)) for p, e in fmpz(n).factor())


def divisors(n):
    result = [1]
    for p, full in factors(n):
        additions = []
        power = p
        while power <= full:
            additions.extend(a*power for a in result)
            power *= p
        result.extend(additions)
    return result


def scan(t, radius):
    hits = {}
    counts = Counter()
    for d in range(-radius, radius+1):
        if d in (0, -1, t-1) or (t == 2 and d > 0):
            continue
        possible_s = divisors(abs((d+1)*(d+1-t)))
        for u in divisors(abs(d)):
            if u < 3:
                continue
            largest_t = isqrt(u**3//27)
            for s in possible_s:
                if gcd(u,s) != 1:
                    continue
                for shift in range((-pow(s,-1,u)) % u, largest_t+1, u):
                    n = 1+s*shift
                    v = n//u
                    if v < 2 or gcd(v,u)>1 or gcd(s,shift)>1:
                        continue
                    a, remainder = divmod(n+d,t)
                    if remainder or not 0 < a < n or a % u:
                        continue
                    counts['assignments'] += 1
                    if any(binomial_carry(n,a,p) for p,w in factors(s)+factors(u)):
                        continue
                    shifted = factors(shift)
                    if any(2*w >= u for p,w in shifted):
                        continue
                    # Assign every prime avoided by this cut to the cut itself.
                    weight = factors(v)
                    if any(not binomial_carry(n,a,p) for p,w in weight+shifted):
                        continue
                    counts['companion tests'] += 1
                    row = find_row(n,v,list(weight),shift,list(shifted))
                    if row is not None:
                        assert sum(row) == n and min(row) > 0
                        assert gcd(n*(n-1), gcd(comb(n,a), comb(n,row[0])*comb(n-row[0],row[1]))) == 1
                        hits[n,a,d] = row
    return counts, hits


def shifted_fractions():
    for denominator, k, p in ((5,3,3),(7,3,3),(7,5,5),(7,6,3)):
        modulus = p*denominator
        residue, = [m for m in range(modulus)
                    if k*m % denominator == 1 and (denominator*m+1) % p == 0]
        graphs = [CarryGraph(base,(k,denominator-k),(0,1)) for base in (2,p,denominator)]
        expected = [8452] if (denominator,k) == (7,5) else []
        accepted, nodes, tested = intersection(graphs,10**1000,modulus,residue)
        assert accepted == expected
        print('Shifted fraction:', k, denominator, accepted, nodes, tested, flush=True)

    n, a = 59165, 42260
    q = n//5
    cut = comb(n,a)
    assert factors(n) == ((5,5),(q,q))
    for i,j,k in ((1,1,3),(1,2,2)):
        assert gcd(n*(n-1),gcd(cut,comb(n,i*q)*comb(n-i*q,j*q))) == 1


def exact_small():
    checks = 0
    for n in range(3,33):
        m = n-1
        coefficients = [comb(n,b)*comb(n-b,c) for b in range(1,n//3+1)
                        for c in range(b,(n-b)//2+1)]
        for a in range(1,n):
            cut = comb(n,a)
            for t in range(2,8):
                for r in range(1,t):
                    if gcd(r,t)>1:
                        continue
                    d = t*a-r*n
                    if d in (0,-r,t-r):
                        continue
                    limit = lcm(abs(d+r), abs(d+r-t))
                    for coefficient in coefficients:
                        weight = gcd(n,gcd(cut,coefficient))
                        shifted = gcd(m,gcd(cut,coefficient))
                        scale_bound = gcd(abs(d),n//weight)
                        shift_bound = gcd(limit,m//shifted)
                        assert 27*m*m <= shifted**2*(weight*scale_bound)**3*shift_bound**2
                        checks += 1
    return checks


if __name__ == '__main__':
    print('Exact retained-divisor checks:', exact_small())
    for t, radius in RADII.items():
        counts, hits = scan(t,radius)
        n, a = (21809,5424) if t == 4 else (5752,568)
        assert set(hits) == {(n,a,t*a-n)} and abs(t*a-n) == radius
        print(t, radius, dict(counts), hits, flush=True)
    shifted_fractions()
