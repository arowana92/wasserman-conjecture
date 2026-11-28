"""Prime lists from a large weight factor and large predecessor powers."""

from itertools import product
from math import comb, prod

from flint import fmpz

from check_heavy_predecessor_crt import factor
from check_direction_pair_count import valuation
from search_endpoint import avoids


def predecessor_rows(n, powers):
    """All positive sorted rows avoiding the specified predecessor bases."""
    shift = prod(w for _,w in powers)
    for positions in product(range(3),repeat=len(powers)):
        residues = [0,0,0]
        for i,(_,w) in zip(positions,powers):
            d = shift//w
            residues[i] += d*pow(d,-1,w)
        a,b,c = [z % shift or shift for z in residues]
        for x in range(a,n//3+1,shift):
            for y in range(x+(b-x)%shift,(n-x)//2+1,shift):
                row = (x,y,n-x-y)
                if row[2] % shift == c % shift and all(avoids(row,p) for p,_ in powers):
                    yield row


def small_cofactor_list(n, report=False):
    q = max(p for p,_ in factor(n))
    r = n
    while not fmpz(r).is_prime():
        r -= 1
    forced = (q,q,n-2*q)
    assert n//q in (3,4) and max(forced) < r
    gap = n-r
    powers = [(p,w) for p,w in factor(n-1) if w >= gap and valuation(n,forced,p)]
    if not powers:
        return None
    # Every near row and the forced row contain every selected predecessor base.
    # A third row avoiding them all has only a small set of CRT residues.
    third = list(predecessor_rows(n,powers))
    primary = {q,r} | {p for p,_ in powers}
    possible = {p for j in range(gap) for p,_ in factor(n-j)}
    extra = sorted(p for p in possible-primary if valuation(n,forced,p))
    constraints = set()
    for k in range(2,gap+1):
        outer = comb(n,k)
        for a in range(1,k//2+1):
            coefficient = outer*comb(k,a)
            eligible = [p for p in extra if coefficient%p == 0]
            for row in third:
                candidates = frozenset(p for p in eligible if valuation(n,row,p))
                if not candidates:
                    return None
                constraints.add(candidates)
    selected = set()
    while constraints:
        p = max(extra,key=lambda p:(sum(p in m for m in constraints),-p))
        selected.add(p)
        constraints = {m for m in constraints if p not in m}
    if report:
        print(n,"gap",gap,"predecessor cofactor",(n-1)//prod(w for _,w in powers),
              "third rows",len(third),"additional primes",sorted(selected))
    return tuple(sorted(primary | selected))
