"""Verify individual near-prime lists by exact coefficients and orthogonal CRT residues."""

import argparse
from collections import Counter
from functools import lru_cache
from math import comb, prod
from pathlib import Path

from flint import fmpz

from check_direction_pair_count import valuation


def certify(n, primes):
    assert primes == sorted(set(primes)) and all(fmpz(p).is_prime() for p in primes)
    if len(primes) == 1:
        p, = primes
        power = 1
        while power*p <= n:
            power *= p
        # Spend the largest digit unit. The rest is zero or one p-power.
        rest = n-power
        if rest == 0:
            return True
        while rest % p == 0:
            rest //= p
        return rest == 1
    *base, near = primes
    assert n < 2*near <= 2*n
    general = any(n*(n-1) % p for p in base)
    atoms = []
    for p in base:
        if n*(n-1) % p:
            atoms.append((p,1,None))
            continue
        m = n if n % p == 0 else n-1
        w = p
        while m % (w*p) == 0:
            w *= p
        atoms.append((p,w,m == n))

    @lru_cache(None)
    def data(mask):
        items = [a for i,a in enumerate(atoms) if (mask >> i) & 1]
        scale = prod(w for _,w,at_weight in items if at_weight)
        powers = [w for _,w,at_weight in items if at_weight is False]
        target = n//scale
        t = prod(powers)
        if target < 3 or any(2*w+1 > target for w in powers) or 27*t*t > target**3:
            return (0,scale,target,t,items,powers)
        return (3**len(powers)*(target//t+1)**2,scale,target,t,items,powers)

    @lru_cache(None)
    def excluded(mask):
        cost,scale,target,t,items,powers = data(mask)
        if cost == 0:
            return True
        if cost > 1000000:
            return False
        # Each part corresponds to a subset of exceptional prime powers.
        # Two normalized parts cannot use the same exceptional power.
        residues = [0]
        for w in powers:
            idempotent = (t//w)*pow(t//w,-1,w)*target
            residues += [(x+idempotent) % t for x in residues]
        for i,ra in enumerate(residues):
            for a in range(ra or t,target//3+1,t):
                for j,rb in enumerate(residues):
                    if i & j:
                        continue
                    for b in range(a+(rb-a)%t,(target-a)//2+1,t):
                        row = (a*scale,b*scale,(target-a-b)*scale)
                        if all(valuation(n,row,p) == 0 for p,_,_ in items):
                            return False
        return True

    # At each gap k the two small parts sum to k. Compute their exact
    # coefficient as a binomial product; take supports independently of Kummer.
    supports = set()
    for k in range(2,n-near+1):
        cut = comb(n,k)
        for a in range(1,k//2+1):
            coefficient = cut*comb(k,a)
            mask = sum(1 << j for j,p in enumerate(base) if coefficient % p == 0)
            supports.add(mask)
    if general:
        # The seven lists with an outside prime also specify q | n with
        # normalized sum 3 or 4. There is only one positive q-avoiding row.
        q = next(p for p in base if n % p == 0 and n//p in (3,4))
        forced = (q,q,n-2*q)
        assert max(forced) < near
        present = sum(1 << j for j,p in enumerate(base) if valuation(n,forced,p))
        q_bit = 1 << base.index(q)
        return all(s & q_bit and excluded(s & present) for s in supports)
    for support in supports:
        subset = support
        while True:
            other = support ^ subset
            if subset <= other:
                sides = sorted((subset,other),key=lambda x:data(x)[0])
                if not any(excluded(x) for x in sides):
                    return False
            if subset == 0:
                break
            subset = (subset-1) & support
    return True


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("certificates")
    args = parser.parse_args()
    stats = Counter()
    seen = set()
    for line in Path(args.certificates).read_text().splitlines():
        n,*primes = map(int,line.split())
        assert n not in seen and certify(n,primes), (n,primes)
        stats[len(primes)] += 1
        seen.add(n)
    print("independently verified",len(seen),"certificate sizes",sorted(stats.items()))
