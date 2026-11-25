"""Classify successor weights p**a*(p**j+d)+1 without ordering the exponents."""

from argparse import ArgumentParser
from bisect import bisect_left, bisect_right
from collections import Counter
from functools import lru_cache
from math import gcd, isqrt, lcm
from typing import NamedTuple

from check_carry_intersections import COFACTORS
from check_close_powers import bases, subproducts
from check_prime_power_cofactors import carry_language
from check_shifted_prefixes import stream_accepts
from check_six_bit_rays import normalized_rows
from check_six_digit_gaps import complement_by_parts
from prime_certificates import primes_through
from search_endpoint import atoms, smallest_factors
from search_even_endpoint import avoids


class Case(NamedTuple):
    Q: int
    u: int
    T: int
    first: tuple
    other: tuple
    p: int
    height: int
    residue: int
    modulus: int


@lru_cache(maxsize=2048)
def powers(base, modulus):
    """A complete period of successive powers of a unit, starting with exponent 0."""
    if gcd(base, modulus) != 1:
        raise ValueError((base, modulus))
    values = [1]
    value = base % modulus
    while value != 1:
        values.append(value)
        value = value*base % modulus
    return values


@lru_cache(maxsize=4096)
def complementary(u, powers_in_T):
    rows = {tuple(sorted(v)) for v in normalized_rows(u, powers_in_T)}
    if rows != complement_by_parts(u, powers_in_T):
        raise ArithmeticError(("complementary rows", u, powers_in_T))
    return sorted(rows)


def first_progression(Q, u, B, factors, x, y):
    # x and y denote the two unshifted slopes in this function.
    residue = -pow(Q, -1, u) % u
    modulus = u
    for r, w in factors:
        if B % r or x % w == y % w == 0:
            continue
        if not ((x % w == 0 and y % r) or (y % w == 0 and x % r)):
            return None
        residue += modulus*((pow(x+y, -1, w)-residue)*pow(modulus, -1, w) % w)
        modulus *= w
    return residue, modulus


def templates(limit, counts):
    maximum = isqrt(limit**5//27)
    if maximum < 4:
        return
    spf = smallest_factors(maximum)
    base_at = {}
    for p in primes_through(maximum):
        height = p
        while height <= maximum:
            base_at[height] = p
            height *= p
    heights = sorted(base_at)
    for Q in range(4, maximum+1):
        factors = atoms(Q, spf)
        if max(w for _, w in factors) > limit:
            continue
        if Q <= 1000 and all(Q % c for c in COFACTORS):
            continue
        available = []
        for H in heights[bisect_left(heights,Q-limit):bisect_right(heights,Q-4)]:
            p = base_at[H]
            cofactor = Q
            while cofactor % p == 0:
                cofactor //= p
            if cofactor > 1000 or (p != 2 and cofactor in COFACTORS):
                available.append((p,H))
        if not available:
            continue
        for T in subproducts(w for _,w in factors if 2*w+1 < limit):
            B = Q//T
            full = tuple(w for r,w in factors if T % r == 0)
            for u in range(2*max(full, default=1)+1,limit):
                if gcd(u,Q) != 1 or 27*T*T > u**3 or B*u >= limit**2:
                    continue
                for y in range(u,limit,u):
                    e = B//gcd(B,y//u)
                    if e > limit or gcd(e,y) != 1:
                        continue
                    step = u*e
                    shifted_small = (e*(Q*pow(e,-1,u) % u)-y) % step
                    for shifted, start in ((False,step),(True,shifted_small or step)):
                        for x in range(start,limit-y+1,step):
                            t = x+y
                            first = (y,Q-t,x) if shifted else (y,x,Q-t)
                            progression = first_progression(Q,u,B,factors,*first[:2])
                            if progression is None:
                                continue
                            residue, modulus = progression
                            counts["first templates"] += 1
                            for p,H in available:
                                if t > Q-H or T % p == 0 or not avoids(Q,first,p):
                                    continue
                                counts["compatible heights"] += 1
                                if gcd(p,modulus) != 1 or residue not in powers(p % modulus,modulus):
                                    continue
                                counts["power orbit"] += 1
                                for other in complementary(u,full):
                                    yield Case(Q,u,T,first,other,p,H,residue,modulus)


def required_bases(c):
    first = set(bases(c.Q//c.T) + bases(c.u)) - {c.p}
    other = set(bases(c.T))
    if c.Q % 2 and c.p % 2 and c.u % 2:
        other.add(2)
    return sorted(first), sorted(other)


def local_parts(c, value, w, first):
    if first:
        x,y,z = c.first
        return x*value, y*value, z*value+1
    v = (c.Q*value+1)*pow(c.u,-1,w) % w
    return tuple(v*f for f in c.other)


def compatible(constraints):
    # Only common prime powers matter for pairwise CRT compatibility.
    projected = []
    for i, (m, allowed) in enumerate(constraints):
        other = lcm(*(n for k,(n,_) in enumerate(constraints) if k != i))
        g = gcd(m,other)
        projected.append((g,{a % g for a in allowed}))
    projected.sort(key=lambda item: len(item[1]))

    def lift(i, stride, residue):
        if i == len(projected):
            return True
        m, allowed = projected[i]
        shared = gcd(stride,m)
        for a in allowed:
            if (a-residue) % shared:
                continue
            step = (a-residue)//shared*pow(stride//shared,-1,m//shared) % (m//shared)
            if lift(i+1,lcm(stride,m),residue+stride*step):
                return True
        return False
    return lift(0,1,0)


def orbit_sieve(c, cap):
    initial = powers(c.p % c.modulus,c.modulus)
    L = len(initial)
    a0 = initial.index(c.residue)
    constraints = [(L,[a0])]
    first, other = required_bases(c)
    for r in first+other:
        w = r
        while w*r <= cap:
            w *= r
        cycle = powers(c.p % w,w)
        stride = gcd(L,len(cycle))
        good = []
        for a in range(a0 % stride,len(cycle),stride):
            value = cycle[a]
            full_u = gcd(c.u,w)
            if 1 < full_u < w and (c.Q*value+1) % (r*full_u) == 0:
                continue
            parts = local_parts(c,value,w,r in first)
            layer = r
            while layer <= w and sum(v % layer for v in parts) < layer:
                layer *= r
            if layer > w:
                good.append(a)
        if not good:
            return False
        constraints.append((len(cycle),good))
    return compatible(constraints)


def affine_sieve(c):
    first, other = required_bases(c)
    for r in first:
        if not stream_accepts(r,c.first,(0,0,1),c.residue,c.modulus,c.Q,c.u):
            return False
    v = (c.Q*c.residue+1)//c.u
    step = c.Q*(c.modulus//c.u)
    for r in other:
        if not carry_language(r,tuple(f*step for f in c.other),tuple(f*v for f in c.other),range(r)):
            return False
    return True


def fixed_divisor(c):
    # For r dividing V, remove its full power from V and look at the first
    # nonzero scalar digit. Try every possible unit digit at this position.
    orbit = powers(c.p % c.modulus,c.modulus)
    L, a0 = len(orbit), orbit.index(c.residue)
    for r in primes_through(c.u):
        if c.u % r == 0 or pow(c.p,L,r) != 1 or (c.Q*pow(c.p,a0,r)+1) % r:
            continue
        if all(sum((v*k) % r for k in c.other) >= r for v in range(1,r)):
            return r
    return None


def layer_candidates(c, limit, base):
    first, other = required_bases(c)
    if base not in first+other:
        raise ValueError("lifting requires a forced base")
    orbit = powers(c.p % c.modulus,c.modulus)
    stride, start = len(orbit), orbit.index(c.residue)
    if pow(c.p,stride,base) != 1:
        raise ValueError("initial period must fix the power at the lifting prime")
    modulus = 1
    levels = []
    while stride <= limit:
        modulus *= base
        branches = 1 if pow(c.p,stride,modulus) == 1 else base
        levels.append((modulus,stride,branches))
        stride *= branches
    # Store the power at final precision. At each split, moving one stride
    # multiplies it by a fixed unit. Depth-first traversal bounds memory.
    levels = [(w,L,b,pow(c.p,L,modulus)) for w,L,b in levels]
    if base in first:
        x,y,z = c.first
        def layer(value,w):
            return x*value % w + y*value % w + (z*value+1) % w < w
    else:
        inverse = pow(c.u,-1,modulus)
        def layer(value,w):
            v = (c.Q*value+1)*inverse % w
            return sum(f*v % w for f in c.other) < w

    def lift(i,a,value):
        if i == len(levels):
            if 1 <= a <= limit:
                yield a
            return
        w,L,branches,multiplier = levels[i]
        for _ in range(branches):
            if a > limit:
                break
            if layer(value,w):
                yield from lift(i+1,a,value)
            a += L
            value = value*multiplier % modulus

    return modulus, lift(0,start,pow(c.p,start,modulus))


def bounded_exponents(c, limit, base=2):
    first, other = required_bases(c)
    _, candidates = layer_candidates(c,limit,base)
    # An independent small binary-layer table filters the leaves cheaply.
    low = None
    if base != 2 and 2 in first+other:
        cycle = powers(c.p % 4096,4096)
        low = [all(sum(k % w for k in local_parts(c,v,w,2 in first)) < w
                   for w in (2**j for j in range(1,13))) for v in cycle]
    for a in candidates:
        if low is not None and not low[a % len(low)]:
            continue
        for r in first+other:
            w = r**100
            entries = local_parts(c,pow(c.p,a,w),w,r in first)
            while w > 1 and sum(k % w for k in entries) < w:
                w //= r
            if w > 1:
                break
        else:
            yield a


def solve(limit, exponent_limit):
    counts = Counter()
    remaining = set()
    for c in templates(limit,counts):
        counts["row pairs"] += 1
        if orbit_sieve(c,256):
            counts["through modulus cap 256"] += 1
            if orbit_sieve(c,4096):
                counts["through modulus cap 4096"] += 1
                if affine_sieve(c):
                    remaining.add(c)
    print(dict(counts), flush=True)
    for c in sorted(remaining):
        print(c,flush=True)
        r = fixed_divisor(c)
        if r:
            counts["fixed divisor"] += 1
            print(f"Excluded at every exponent by the fixed divisor {r}.",flush=True)
            continue
        base = 7 if c.u == 7 else 2
        survivors = list(bounded_exponents(c,exponent_limit,base))
        if survivors:
            raise ArithmeticError(("unresolved in bounded interval",c,survivors[:5]))
    print(f"{len(remaining)-counts['fixed divisor']} remaining templates, "
          f"none possible at exponents up to {exponent_limit}.")


if __name__ == "__main__":
    cli = ArgumentParser(description=__doc__)
    cli.add_argument("max_offset",type=int,nargs="?",default=256)
    cli.add_argument("--exponent-limit",type=int,default=10**24)
    args = cli.parse_args()
    if not 1 <= args.max_offset <= 256 or args.exponent_limit < 1:
        cli.error("use 1 through 256 for the offset and a positive exponent limit")
    solve(args.max_offset,args.exponent_limit)
