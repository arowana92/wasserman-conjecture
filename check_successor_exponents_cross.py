"""Independently check selected unrestricted-exponent row templates and the exponent sieves."""

from collections import Counter
from functools import lru_cache
from itertools import product
from math import comb, lcm

from carry_blocks import CarryGraph

from check_successor_exponents import Case, bounded_exponents, compatible, fixed_divisor, layer_candidates, orbit_sieve, templates
from search_endpoint import avoids
from search_even_endpoint import prime_powers


def check_row_pairs():
    targets = {(364,11,4,283),(1651,4,1,1459),(5126,7,2,4889)}
    forward = {target: set() for target in targets}
    for c in templates(256,Counter()):
        key = (c.Q,c.u,c.T,c.p)
        if key in targets and c.height == c.p:
            forward[key].add((c.first,c.other))
    for Q,u,T,p in sorted(targets):
        wlist = list(prime_powers(Q//T).values())
        small = Q-p
        other = []
        for f in range(1,u//3+1):
            for g in range(f,(u-f)//2+1):
                row = (f,g,u-f-g)
                if all(sorted(k % w for k in row) == [0,0,u % w]
                       for w in prime_powers(T).values()):
                    other.append(row)
        m = u*(Q//T)
        values = set()
        value = p % m
        while value not in values:
            values.add(value)
            value = value*p % m
        allowed = {v for v in values if (Q*v+1) % u == 0}

        @lru_cache(None)
        def residues(w,a,b):
            z = (Q-a-b) % w
            return {v for v in allowed if sorted((v*a % w,v*b % w,(v*z+1) % w)) == [0,0,1]}

        direct = set()
        for y in range(u,small,u):
            for x in range(1,small-y+1):
                z = Q-x-y
                for first in ((y,x,z),(y,z,x)):
                    a,b,c = first
                    if b % u or not avoids(first,p):
                        continue
                    possible = allowed.copy()
                    for w in wlist:
                        possible.intersection_update(residues(w,a % w,b % w))
                    if possible:
                        direct.update((first,row) for row in other)
        if not direct or direct != forward[Q,u,T,p]:
            raise ArithmeticError(("template comparison",Q,direct ^ forward[Q,u,T,p]))
        print(f"Direct residue enumeration agrees on {len(direct)} row pairs at Q={Q}, p={p}.",flush=True)


def check_sieves():
    # The small examples exercise local filters in isolation.
    for m,n,k in product(range(2,9),repeat=3):
        sets = [(v,[j for j in range(v) if (j*j+i) % 5 < 2]) for i,v in enumerate((m,n,k))]
        direct = any(all(t % v in vals for v,vals in sets) for t in range(lcm(m,n,k)))
        if direct != compatible(sets):
            raise ArithmeticError(("CRT",sets))
    hits = 0
    for x,y in ((4,16),(4,4),(8,12),(8,16)):
        c = Case(27,4,1,(x,y,27-x-y),(1,1,2),3,3,1,4)
        expected = []
        for a in range(1,19):
            P = 3**a
            parts = (P*x,P*y,P*(27-x-y)+1)
            if P % 4 == 1 and avoids(parts,2):
                expected.append(a)
        if list(bounded_exponents(c,18)) != expected:
            raise ArithmeticError(("binary lift",c,expected))
        if expected:
            hits += len(expected)
        # Period of powers of 3 modulo 16 is four. Require the full 2-part of N to be 4.
        good = False
        for a in range(4):
            P = 3**a
            row = (P*x,P*y,P*(27-x-y)+1)
            good |= P % 4 == 1 and (27*P+1) % 8 == 4 and all(
                sum(row)//w == sum(v//w for v in row) for w in (2,4,8,16))
        if orbit_sieve(c,16) != good:
            raise ArithmeticError(("modular orbit",c))
    if not hits or comb(244,36)*comb(208,144) % 2 != 1:
        raise ArithmeticError("lost the weight-244 positive control")
    print(f"343 generalized CRT comparisons and {hits} positive exponent controls agree.")


def compare_lifts():
    for c,r,a0,L in (
        (Case(1651,4,1,(104,1524,23),(1,1,2),1459,1459,2297,6604),2,86,126),
        (Case(5126,7,2,(231,4893,2),(1,2,4),4877,4877,8970,17941),7,587,870),
    ):
        limit = 10**6
        modulus, actual = layer_candidates(c,limit,r)
        expected = set()
        for a in range(a0,limit+1,L):
            P = pow(c.p,a,modulus)
            x,y,z = c.first
            row = (P*x,P*y,P*z+1)
            w = r
            while w <= modulus and sum(row)//w == sum(k//w for k in row):
                w *= r
            if w > modulus:
                expected.add(a)
        if not expected or expected != set(actual):
            raise ArithmeticError(("lifting classes",r))
        print(f"Direct exponents agree on {len(expected)} survivors of the base-{r} lift.")
    c = Case(5126,7,2,(231,4893,2),(1,2,4),4889,4889,8970,17941)
    if fixed_divisor(c) != 5:
        raise ArithmeticError("the common prime 5 was missed")


def check_path_counts():
    for r,slopes,n,num,den in ((2,(104,1524,23),1024,1,5),(7,(231,4893,2),1024,5,16)):
        graph = CarryGraph(r,slopes,(0,0,1))
        paths = graph.block_bound(n)
        if paths**den >= r**(num*n):
            raise ArithmeticError(("path bound",r,paths))
        print(f"{len(graph.states)} states, base {r}, certified exponent at most {num}/{den}.")


if __name__ == "__main__":
    check_sieves()
    compare_lifts()
    check_row_pairs()
    check_path_counts()
