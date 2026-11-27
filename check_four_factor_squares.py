"""Exhaust local restrictions on the two nonendpoint rows at 1 + D t^2."""

from functools import lru_cache
from math import gcd

from check_three_factors import powers


def digits_ok(row, p, depth):
    for _ in range(depth):
        if sum(a % p for a in row) >= p:
            return False
        row = tuple(a//p for a in row)
    return True


@lru_cache(None)
def squares(modulus):
    return {t*t % modulus for t in range(modulus)}


def local(D, x, U, f, prime, assigned_power):
    depth = 3
    modulus = prime**depth
    while modulus <= assigned_power:
        modulus *= prime
        depth += 1
    for W in squares(modulus):
        n = 1+D*W
        cut = (x*W, (D-x)*W+1)
        if assigned_power:
            if n % assigned_power or n % (prime*assigned_power) == 0:
                continue
            if digits_ok(cut,prime,depth):
                return True
        else:
            # D is squarefree. If the base is also a base of t, its full
            # predecessor power contributes one extra copy of prime to x.
            if W % prime == 0 and x % prime:
                continue
            V = n*pow(U,-1,modulus) % modulus
            if digits_ok(cut,prime,depth) or digits_ok(tuple(V*a for a in f),prime,depth):
                return True
    return False


def candidates(D):
    result = []
    for U in range(3,D):
        if gcd(U,D) != 1:
            continue
        primes_U = powers(U)
        for x in range(U,D,U):
            if not all(local(D,x,U,(),p,w) for p,w in primes_U.items()):
                continue
            for a in range(1,U//3+1):
                for b in range(a,(U-a)//2+1):
                    f = a,b,U-a-b
                    if all(local(D,x,U,f,p,0) for p in powers(D)):
                        result.append((x,U,f))
    return result


def independent(D):
    # Use the exact characterization of squares modulo a prime power and
    # independent factorial floors, traversing normalized triples first.
    @lru_cache(None)
    def residues(p, k):
        modulus = p**k
        out = [0]
        for W in range(1,modulus):
            z, v = W, 0
            while z%p == 0:
                z //= p
                v += 1
            if v%2:
                continue
            if p == 2:
                ok = z % min(8,2**(k-v)) == 1
            else:
                ok = pow(z,(p-1)//2,p) == 1
            if ok:
                out.append(W)
        return tuple(out)

    @lru_cache(None)
    def permissible(p, x, U, shape):
        power = 1
        while U % (p*power) == 0:
            power *= p
        k = 3
        while p**k <= power:
            k += 1
        levels = tuple(p**j for j in range(1,k+1))
        mod = levels[-1]
        for W in residues(p,k):
            N = D*W+1
            part = x*W
            cut_cost = sum(N//w-part//w-(N-part)//w for w in levels)
            if power > 1:
                if N%power == 0 and N%(p*power) and cut_cost == 0:
                    return True
                continue
            if W%p == 0 and x%p:
                continue
            if cut_cost == 0:
                return True
            V = N*pow(U,-1,mod) % mod
            if all(U*V//w == sum(a*V//w for a in shape) for w in levels):
                return True
        return False

    out = set()
    for a in range(1,D):
        for b in range(a,D):
            for c in range(b,D-a-b):
                U = a+b+c
                if gcd(U,D) != 1:
                    continue
                shape = (a,b,c)
                for x in range(U,D,U):
                    if (all(permissible(p,x,U,()) for p in powers(U))
                        and all(permissible(p,x,U,shape) for p in powers(D))):
                        out.add((x,U,shape))
    return out


def fourth_powers():
    # This single modulus tests all the unit residue classes used in the two
    # hand proofs. Primes dividing the modulus are checked separately.
    L = 16*3*5*7*11*17
    units = {pow(z,4,L) for z in range(1,L,2) if gcd(z,L)==1}
    def any_candidate(D,R,b,r=1):
        for x in range(3,D):
            if x*R%16 & ~(D*R%16):
                continue
            if D%r == 0 and x%r:
                continue
            for U in range(3,x+1,2):
                if x%U or gcd(U,D*r)>1:
                    continue
                W = (1<<b)*R
                if (1+D*W)%U:
                    continue
                # Each rejected static base must be taken by the companion.
                rejected = [p for p in powers(D) if (x*W)%p not in (0,1)]
                if any(2*p+1>U for p in rejected):
                    continue
                if len(rejected)==2 and D>=U:
                    continue
                return True
        return False
    # For unit and nonunit R the low binary test always leaves only U=5,17,
    # or U=17 respectively, so the order of 2 at all required odd moduli divides 120.
    tests = 0
    for D in (21,33):
        for R in units:
            for b in range(1,121):
                assert not any_candidate(D,R,b)
                tests += 1
        for r in (3,5,7,11,17):
            R, seen = pow(r,4,L), set()
            while R not in seen:
                seen.add(R)
                for b in range(1,121):
                    assert not any_candidate(D,R,b,r)
                    tests += 1
                R = R*pow(r,4,L)%L
        R, seen = 17%L, set()
        while R not in seen:
            seen.add(R)
            assert all(not any_candidate(D,R,b,17) for b in range(1,121))
            R = R*17%L
    return tests


if __name__ == "__main__":
    for D in (21,33):
        rows = candidates(D)
        assert set(rows) == independent(D)
        for x,U,shape in rows:
            assert U%2 and (x & ~D & 7)
            V = (1+D)*pow(U,-1,8) % 8
            assert not digits_ok(tuple(V*a for a in shape),2,3)
        print(D, rows)
    print("fourth-power residue checks", fourth_powers())
