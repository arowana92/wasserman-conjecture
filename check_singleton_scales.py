"""Classify an endpoint row with all assigned predecessor powers exceptional in one slot."""

from argparse import ArgumentParser
from collections import Counter, deque
from functools import lru_cache
from math import gcd, isqrt

from check_three_factors import powers
from search_endpoint import avoids


@lru_cache(maxsize=32768)
def atoms(n):
    return tuple(powers(n).items())


def companions(limit):
    """Choose the three coprime exceptional-position factors in increasing order."""
    for d in range(3, isqrt(limit//3)+1):
        for e in range(d+1, limit):
            if d*e + (d+e)*e >= limit:
                break
            if gcd(d,e) != 1:
                continue
            for h in range(e+1, (limit-d*e)//(d+e)+1):
                if gcd(h,d*e) != 1:
                    continue
                T = d*e*h
                cap = min(limit,(T-1)//2)
                for x in range(e*h,cap-d*h-d*e+1,e*h):
                    for y in range(d*h,cap-x-d*e+1,d*h):
                        for z in range(d*e,cap-x-y+1,d*e):
                            U = x+y+z
                            if gcd(U,T) == 1:
                                yield U,T,tuple(sorted((x,y,z)))


@lru_cache(maxsize=150000)
def accepts_prefix(q, slopes, offsets, step, residue, modulus, forbidden):
    """Decide existence at the prescribed low digits by exact graph reachability."""
    states = set()
    for prefix in range(residue, modulus, step):
        if prefix == forbidden:
            continue
        row = tuple(a*prefix+b for a,b in zip(slopes,offsets))
        layer = 1
        while layer < modulus:
            layer *= q
            if sum(x % layer for x in row) >= layer:
                break
        else:
            states.add(tuple(x//modulus for x in row))
    queue = deque(states)
    while queue:
        state = queue.popleft()
        if avoids(state,q):
            return True
        for digit in range(q):
            row = tuple(a*digit+x for a,x in zip(slopes,state))
            if sum(x % q for x in row) >= q:
                continue
            dest = tuple(x//q for x in row)
            if dest not in states:
                states.add(dest)
                queue.append(dest)
    return False


def assigned(q, w, slopes, offsets, coefficient, sign):
    residue = sign*pow(coefficient,-1,w) % w
    forbidden = sign*pow(coefficient,-1,q*w) % (q*w)
    return accepts_prefix(q,slopes,offsets,w,residue,q*w,forbidden)


ZERO = (0,0,0)
SHIFT = (1,0,0)


def surviving_shapes(limit, counts):
    for U,T,y in companions(limit):
        counts['companion systems'] += 1
        if not all(assigned(q,w,y,ZERO,U,1) for q,w in atoms(T)):
            continue
        counts['companion systems passing local tests'] += 1
        for b in range(2,(T-1)//U+1):
            k = T-U*b
            for f in range(1,b//2+1):
                g = b-f
                slopes = (k,U*f,U*g)
                if not all(assigned(q,w,slopes,SHIFT,T,-1) for q,w in atoms(U)):
                    continue
                if T%2 and U%2:
                    # At even S the new units bit is separate in the shifted row.
                    # At even V the ordinary row can be divided by its power of 2.
                    if not any(accepts_prefix(2,a,ZERO,2,1,2,-1) for a in (slopes,y)):
                        continue
                counts['surviving shapes'] += 1
                yield U,T,k,f,g,y


def pure_power_test(shape, counts):
    U,T,k,f,g,y = shape
    if len(atoms(U)) != 1:
        return
    p = atoms(U)[0][0]
    counts['prime-power scales'] += 1
    W = U
    while W < p**16:
        if (W-1) % k == 0:
            S = (W-1)//k
            V = (T*S+1)//U
            if (gcd(S,T) == gcd(V,U) == 1 and avoids((W,U*f*S,U*g*S),p)
                    and all(avoids(tuple(V*x for x in y),q) for q,_ in atoms(T))):
                raise ArithmeticError(('small exponent', shape, W))
        W *= p
    if T%2 == 0 and k%2 and p%2:
        # S is even whenever kS+1 is an odd prime power.
        counts['parity exclusions'] += 1
        return
    modulus = p**16
    S0 = -pow(k,-1,modulus) % modulus
    x,z = f*S0 % modulus, g*S0 % modulus
    assert not avoids((x,z,0),p), ('undecided large exponent',shape)
    counts['low digit exclusions'] += 1


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--limit', type=int, default=5000)
    args = parser.parse_args()
    counts = Counter()
    for shape in surviving_shapes(args.limit,counts):
        pure_power_test(shape,counts)
        print('SHAPE',*shape)
    print(dict(counts))


if __name__ == '__main__':
    main()
