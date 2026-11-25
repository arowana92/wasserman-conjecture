"""Check cofactor <=1000 for all small prime bases, with an exponent bound on surviving rays."""

from argparse import ArgumentParser
from math import gcd

from carry_blocks import CarryGraph
from check_carry_intersections import dimension_percent
from check_direction_prime_bases import (
    SPF, every_unit_carries, local_targets, no_low_carry, odd_shapes, orbit_classes,
)
from check_successor_exponents import Case, compatible, layer_candidates, powers
from search_endpoint import atoms, avoids, smallest_factors


def local_conditions(m, shape, residue, modulus, cap):
    u = shape[0]
    m_factors = dict(atoms(m,SPF))
    u_factors = dict(atoms(u,SPF))
    modulus_factors = dict(atoms(modulus,SPF))
    primes = {2} | m_factors.keys() | u_factors.keys()
    constraints = []
    for p in range(2,u+1):
        if p not in primes and SPF[p] == p and every_unit_carries(shape[4:],p):
            constraints.append((p,p,set(range(1,p))-{-pow(m,-1,p) % p}))
    for p in sorted(primes):
        w = modulus_factors.get(p,1)
        q = p
        while q*p <= cap:
            q *= p
        q = max(q,w,p*u_factors.get(p,1))
        constraints.append((p,q,local_targets(m,shape,p,q,residue,w)))
    return constraints


def possible_exponents(prime, first, period, conditions):
    classes = [(period,[first])]
    for p,q,targets in conditions:
        if p == prime:
            continue
        cycle = powers(prime % q,q)
        stride = gcd(period,len(cycle))
        allowed = [i for i in range(first % stride,len(cycle),stride) if cycle[i] in targets]
        if not allowed:
            return False
        classes.append((len(cycle),allowed))
    return compatible(classes)


def exponent_check(m,shape,prime,residue,modulus,forced,limit):
    u, a, b, c, *other = shape
    # The shared lifting routine reads only the chosen first-row base.
    # Include every unforced predecessor power in the upper bound on T.
    upper_t = 1
    for p,w in atoms(m,SPF):
        if p not in forced:
            upper_t *= w
    case = Case(m,u,upper_t,(a,b,c),tuple(other),prime,0,residue,modulus)
    def rank(p):
        graph = CarryGraph(p,(a,b,c),(0,0,1))
        return dimension_percent(graph),p
    base = min((p for p in forced if modulus % p == 0),key=rank)
    _, candidates = layer_candidates(case,limit,base)
    low_tables = []
    for p in forced:
        if p == base:
            continue
        q = p
        while q*p <= 65536:
            q *= p
        cycle = powers(prime % q,q)
        low_tables.append(bytes(no_low_carry((a*v,b*v,c*v+1),p,q) for v in cycle))
    leaves = 0
    for exponent in candidates:
        leaves += 1
        if any(not table[exponent % len(table)] for table in low_tables):
            continue
        for p in sorted(forced,key=lambda q:q == base):
            precision = p**60
            value = pow(prime,exponent,precision)
            if not no_low_carry((a*value,b*value,c*value+1),p,precision):
                break
        else:
            raise ArithmeticError(("exponent requires more tests",m,prime,exponent))
    return base,leaves


def verify(maximum,limit):
    spf = smallest_factors(maximum)
    primes = [p for p in range(3,maximum+1) if spf[p] == p]
    total = leaves = rows = 0
    survivors = []
    for m,shape,residue,modulus,forced in odd_shapes():
        allowed = {z:(k,order) for z,k,order in orbit_classes(modulus,residue)}
        bases = [p for p in primes if p % modulus in allowed and m % p and avoids(shape[1:4],p)]
        total += len(bases)
        for cap in (256,4096,65536):
            if not bases:
                break
            conditions = local_conditions(m,shape,residue,modulus,cap)
            bases = [p for p in bases if possible_exponents(p,*allowed[p % modulus],conditions)]
        for prime in bases:
            survivors.append((m,prime))
            base,count = exponent_check(m,shape,prime,residue,modulus,forced,limit)
            leaves += count
            print(f"d={m}, prime={prime}, lifting base {base}, {count} leaves",flush=True)
        rows += 1
    print(f"{rows} progressions, {total} initial base/shape pairs, {len(survivors)} "
          f"unbounded pairs, {leaves} leaves through exponent {limit}.")
    print("Possible pairs at larger exponents:",survivors)


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--prime-limit",type=int,default=10000)
    parser.add_argument("--exponent-limit",type=int,default=10**24)
    args = parser.parse_args()
    if args.prime_limit < 3 or args.exponent_limit < 1:
        parser.error("use a prime bound of at least three and a positive exponent bound")
    verify(args.prime_limit,args.exponent_limit)
