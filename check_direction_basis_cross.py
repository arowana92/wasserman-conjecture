"""Verify the prime-base density by discrete logarithms and direct factorial layers."""

from fractions import Fraction
from functools import lru_cache
from math import gcd, prod

from check_direction_prime_bases import SPF, odd_shapes
from check_prime_power_cofactors import carry_language
from search_endpoint import atoms, smallest_factors


@lru_cache(maxsize=32768)
def finishing(p, slopes, quotients):
    return carry_language(p,slopes,quotients,range(p))


def prefix_extends(t, q, p, slopes, shifted):
    parts = tuple(a*t+int(shifted and i == 2) for i,a in enumerate(slopes))
    layer = p
    while layer <= q:
        if sum(parts)//layer != sum(x//layer for x in parts):
            return False
        layer *= p
    return finishing(p,slopes,tuple(x//q for x in parts))


def conditions(m,shape,r,M):
    u,a,b,c,*others = shape
    others = tuple(others)
    primes = {2} | {v for v,_ in atoms(m,SPF)+atoms(u,SPF)}
    for p in range(2,max(primes | {u})+1):
        if SPF[p] != p:
            continue
        block_c = (all(sum(x*i % p for x in others) >= p for i in range(1,p))
                   or any((x+y)//gcd(x,y) % p == 0 for i,x in enumerate(others)
                          for y in others[i+1:]))
        if p not in primes:
            # A pair binomial also forces divisibility in the large calculation.
            # Every extra odd prime in these six C shapes passes the unit-digit test.
            if block_c:
                yield p,1,p,{i for i in range(1,p) if (m*i+1) % p}
            continue
        w = dict(atoms(M,SPF)).get(p,1)
        q = w*p*p
        du = gcd(u,q)
        mp = dict(atoms(m,SPF)).get(p,1)
        allowed = set()
        impossible_a = (a+b)//gcd(a,b) % p == 0
        possible_a = mp == 1 or (
            max(gcd(a,mp),gcd(b,mp)) == mp and min(gcd(a,p),gcd(b,p)) == 1)
        possible_c = mp == 1 or sum(y % mp == 0 for y in others) >= 2
        for t in range(r % w,q,w):
            if t % p == 0:
                continue
            # Solve U*V = m*t+1 modulo q by enumerating its reduced congruence.
            v = (m*t+1)//du*pow(u//du,-1,q//du) % (q//du)
            if u % p == 0 and v % p == 0:
                continue
            ok_a = possible_a and not impossible_a and prefix_extends(t,q,p,(a,b,c),True)
            ok_c = possible_c and not block_c and prefix_extends(v,q//du,p,others,False)
            if p == 2 and m % 2 and u % 2:
                if not ok_c:
                    continue
            if (ok_a if u % p == 0 else (ok_a or ok_c)):
                allowed.add(t)
        yield p,w,q,allowed


@lru_cache(maxsize=128)
def log_table(p,q):
    phi = q-q//p
    factors = [v for v,_ in atoms(phi,smallest_factors(phi))]
    root = next(x for x in range(2,q) if x % p and all(pow(x,phi//v,q) != 1 for v in factors))
    logs = [-1]*q
    t = 1
    for i in range(phi):
        logs[t] = i
        t = t*root % q
    assert t == 1 and logs.count(-1) == q-phi
    return logs


def fraction_counter(p,w,q,allowed):
    phi = q-q//p
    if p != 2:
        logs = log_table(p,q)
        good_logs = {logs[t] for t in allowed}
    projections = {}

    @lru_cache(None)
    def count(z,k,order):
        total = matches = 0
        for y in range(z,q,w):
            if y % p == 0:
                continue
            total += 1
            if p == 2:
                # The unit group at 2 is generally not cyclic.
                hits = any(pow(y,k+j*order,q) in allowed for j in range(phi))
            else:
                g = gcd(order*logs[y],phi)
                if g not in projections:
                    projections[g] = {a % g for a in good_logs}
                hits = k*logs[y] % g in projections[g]
            matches += hits
        return Fraction(matches,total)
    return count


def density(m,shape,r,M):
    phi = prod(w-w//p for p,w in atoms(M,SPF))
    factors = atoms(phi,SPF)
    local = [(w,q-q//p,fraction_counter(p,w,q,values))
             for p,w,q,values in conditions(m,shape,r,M)]
    mass = Fraction()
    for z in range(1,M):
        if gcd(z,M) != 1:
            continue
        order = phi
        for p,_ in factors:
            while order % p == 0 and pow(z,order//p,M) == 1:
                order //= p
        # Baby steps and giant steps give k with z^k = r mod M.
        width = 1
        while width*width < order:
            width += 1
        baby = {pow(z,i,M):i for i in range(width)}
        shift = pow(z,-width,M)
        value = r
        for j in range((order+width-1)//width):
            if value in baby:
                k = (j*width+baby[value]) % order
                assert pow(z,k,M) == r
                mass += prod(count(z % w,k % n,order % n) for w,n,count in local)
                break
            value = value*shift % M
    return mass/phi


if __name__ == "__main__":
    upper = sum(density(m,shape,r,M) for m,shape,r,M,_ in odd_shapes())
    print("Independent density sum:",upper,flush=True)
    assert upper == Fraction(14701541568851133613999,530916553722651093597600)
