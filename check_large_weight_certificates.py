"""Find and verify lists of common primes for triples of rows at large weights."""

import argparse
from collections import Counter
from functools import lru_cache
from itertools import combinations
from math import prod
from pathlib import Path

from flint import fmpz

from check_heavy_predecessor_crt import factor
from check_direction_pair_count import valuation
from check_small_cofactor_certificates import small_cofactor_list
from search_endpoint import find_row, splits


def near_masks(n, r, atoms):
    """Inclusion-minimal prime supports of positive rows avoiding r > n/2."""
    gap = n-r
    masks = set()
    for a in range(1,gap//2+1):
        for b in range(a,gap-a+1):
            row = (a,b,n-a-b)
            mask = frozenset(atom for atom in atoms if valuation(n,row,atom[0]) > 0)
            masks.add(mask)
    return [mask for mask in masks if not any(m < mask for m in masks)]


def universal_prime(n):
    # A universal prime divides the endpoint. A carry-free positive split
    # is possible exactly when there are at least three base-p digit units.
    for p,_ in factor(n)+factor(n-1):
        digit_sum = 0
        value = n
        while value and digit_sum < 3:
            value, digit = divmod(value,p)
            digit_sum += digit
        if digit_sum <= 2:
            return p
    return None


def search(n, work_limit, exhaustive):
    single = universal_prime(n)
    if single is not None:
        return (single,)
    r = n
    while not fmpz(r).is_prime():
        r -= 1
    weight = factor(n)
    atoms = tuple(a for a in weight+factor(n-1) if a[0] != r)
    masks = near_masks(n,r,atoms)

    @lru_cache(None)
    def cost(ua,sa):
        scale = prod(w for _,w in ua)
        shift = prod(w for _,w in sa)
        v = n//scale
        if v < 3 or 27*shift**2 > v**3 or any(2*w+1 > v for _,w in sa):
            return 0
        return 3**len(sa)*(v//shift+1)**2

    @lru_cache(None)
    def ruled_out(ua,sa):
        count = cost(ua,sa)
        if count == 0:
            return True
        if count > work_limit:
            return False
        return find_row(n,prod(w for _,w in ua),list(ua),
                        prod(w for _,w in sa),list(sa)) is None

    @lru_cache(None)
    def certifies_pair(support):
        wa = [a for a in atoms if a in support and a in weight]
        ma = [a for a in atoms if a in support and a not in weight]
        for _,ua,va in splits(wa):
            for _,sa,ta in splits(ma):
                sides = [(tuple(ua),tuple(sa)),(tuple(va),tuple(ta))]
                if not any(ruled_out(*side) for side in sorted(sides,key=lambda x:cost(*x))):
                    return False
        return True

    def certifies(choice):
        return all(certifies_pair(mask.intersection(choice)) for mask in masks)

    if not certifies(atoms):
        largest = max(p for p,_ in weight)
        if n//largest in (3,4):
            return small_cofactor_list(n)
        return None
    if exhaustive:
        for k in range(len(atoms)+1):
            for choice in combinations(atoms,k):
                if certifies(choice):
                    return tuple(sorted([r]+[p for p,_ in choice]))
    choice = list(atoms)
    for atom in sorted(atoms,key=lambda x:x[1]):
        trial = [x for x in choice if x != atom]
        if certifies(trial):
            choice = trial
    return tuple(sorted([r]+[p for p,_ in choice]))


def run(path, exhaustive, limit, output):
    weights = {sum(line[3:6])*line[6] for line in
               (tuple(map(int,s.split())) for s in Path(path).read_text().splitlines())}
    hist = Counter()
    with Path(output).open('w') as stream:
        for n in sorted(weights):
            certificate = search(n,limit,exhaustive)
            if certificate is None:
                hist[0] += 1
            else:
                hist[len(certificate)] += 1
                print(n,*certificate,file=stream)
    print("weights",len(weights),"certificate sizes (zero is inconclusive)",sorted(hist.items()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("rows")
    parser.add_argument("output")
    parser.add_argument("--work-limit",type=int,default=1000000)
    parser.add_argument("--exhaustive",action="store_true")
    args = parser.parse_args()
    run(args.rows,args.exhaustive,args.work_limit,args.output)
