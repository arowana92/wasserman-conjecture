"""Regenerate all singleton-position systems by factoring T before choosing the rows."""

from collections import Counter
from functools import lru_cache
from itertools import product
from math import gcd, isqrt, prod

from check_prime_power_cofactors import carry_language
from check_singleton_scales import companions, surviving_shapes
from check_three_factors import powers, valuation
from search_endpoint import atoms, smallest_factors


def denominators(limit, spf):
    for T in range(7, isqrt(limit**3//27)+1):
        ws = [w for _,w in atoms(T,spf)]
        for rest in product(range(3), repeat=len(ws)-1):
            if 1 not in rest or 2 not in rest or rest.index(1) > rest.index(2):
                continue
            groups = (0, *rest)
            d = [prod(w for w,j in zip(ws,groups) if j == i) for i in range(3)]
            if min(d) < 3:
                continue
            m,n,r = (T//v for v in d)
            cap = min(limit,(T-1)//2)
            for x in range(m,cap-n-r+1,m):
                for y in range(n,cap-x-r+1,n):
                    for z in range(r,cap-x-y+1,r):
                        U = x+y+z
                        if gcd(U,T) == 1:
                            yield U,T,tuple(sorted((x,y,z)))


@lru_cache(maxsize=150000)
def check_base(q,w,coefficient,sign,slopes,offsets):
    r = (sign*pow(coefficient,-1,w))%w
    row = tuple(a*r+b for a,b in zip(slopes,offsets))
    layer = q
    while layer <= w:
        if sum(row)//layer > sum(x//layer for x in row):
            return False
        layer *= q
    first = tuple(d for d in range(q) if (coefficient*(r+w*d)-sign)%(q*w))
    return carry_language(q,slopes,tuple(x//w for x in row),first)


def pure(U,T,k,f,g,y,p):
    a,W = 1,p
    while W < U:
        a += 1
        W *= p
    for i in range(a,16):
        W = p**i
        if (W-1)%k:
            continue
        S = (W-1)//k
        weight = T*S+1
        assert weight%U == 0
        V = weight//U
        if gcd(S,T)*gcd(U,V) > 1:
            continue
        A = (W,U*f*S,U*g*S)
        C = tuple(V*x for x in y)
        assert valuation(weight,A,p) or any(valuation(weight,C,q) for q in powers(T))
    if T%2 == 0 and k%2 and p%2:
        return 0
    for j in range(1,17):
        w = p**j
        s = -pow(k,-1,w)
        if (f*s)%w + (g*s)%w >= w:
            return j
    raise AssertionError((U,T,k,f,g,y))


def main():
    limit = 5000
    spf = smallest_factors(isqrt(limit**3//27))
    pending = set(companions(limit))
    initial_count = len(pending)
    found = set()
    max_layer = 0
    for U,T,y in denominators(limit,spf):
        pending.remove((U,T,y))
        if not all(check_base(q,w,U,1,y,(0,0,0)) for q,w in atoms(T,spf)):
            continue
        for f in range(1,T//(2*U)+1):
            for g in range(f,(T-1)//U-f+1):
                k = T-U*(f+g)
                slopes = (k,U*f,U*g)
                if not all(check_base(q,w,T,-1,slopes,(1,0,0)) for q,w in atoms(U,spf)):
                    continue
                if U*T%2 and not any(carry_language(2,z,(0,0,0),[1]) for z in (slopes,y)):
                    continue
                shape = (U,T,k,f,g,y)
                found.add(shape)
                if len(atoms(U,spf)) == 1:
                    max_layer = max(max_layer,pure(*shape,atoms(U,spf)[0][0]))
    assert not pending and found == set(surviving_shapes(limit,Counter()))
    print(f'{initial_count} companion systems, {len(found)} identical shapes; '
          f'smallest remaining scale {min(z[0] for z in found)}; '
          f'largest prime-power carry layer needed {max_layer}.')


if __name__ == '__main__':
    main()
