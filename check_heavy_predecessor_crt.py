"""Verify primitive CRT rows, their endpoint gcds and the full conjecture at their weights."""

import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from math import gcd, prod
from pathlib import Path

from flint import fmpz

from check_direction_pair_count import valuation
from search_endpoint import find_row, splits


@lru_cache(maxsize=50000)
def factor(n):
    result = tuple((int(p), int(p)**int(e)) for p,e in fmpz(n).factor())
    assert prod(w for _,w in result) == n
    assert all(fmpz(p).is_prime() for p,_ in result)
    return result


def direct_box(power_bound, scale_bound):
    """Use the normalized weight first, then independently recover CRT idempotents."""
    atoms = [(m, factor(m)[0][0]) for m in range(2,power_bound+1) if len(factor(m)) == 1]
    hits = set()
    for triple in combinations(atoms, 3):
        powers, bases = zip(*triple)
        if len(set(bases)) != 3:
            continue
        t = prod(powers)
        idempotents = [(t//d)*pow(t//d,-1,d) for d in powers]
        for v in range(2,scale_bound+1):
            if gcd(t,v) != 1:
                continue
            u = pow(v,-1,t)
            y = tuple(u*e % t for e in idempotents)
            if sum(y) != u or gcd(*y) != 1:
                continue
            n = u*v
            row = tuple(x*v for x in y)
            if all(valuation(n,row,p) == 0 for p in bases):
                hits.add((*powers,*y,v))
    return hits


def endpoint_gcd(n, row):
    result = 1
    avoided = 1
    for m in (n,n-1):
        for p,w in factor(m):
            exponent = valuation(n,row,p)
            result *= gcd(w, p**exponent)
            if not exponent:
                avoided *= w
    return result, avoided


def check_full_weight(n, stats):
    weight, pred = factor(n), factor(n-1)
    # One side of an endpoint assignment receives S >= normalized sum V.
    heavy = [(s,sa,ta) for s,sa,ta in splits(pred) if len(sa) >= 3]
    for u,ua,va in splits(weight):
        v = n//u
        if min(u,v) < 3:
            continue
        for s,sa,ta in heavy:
            if s < v or 27*s*s > v**3 or any(2*w+1 > v for _,w in sa):
                continue
            stats['heavy assignments'] += 1
            first = find_row(n,u,ua,s,sa)
            if first is None:
                continue
            stats['admissible first assignments'] += 1
            second = find_row(n,v,va,(n-1)//s,ta)
            if second is not None:
                raise ArithmeticError((n,first,second))


def direct_companion(n,row):
    # A companion avoiding the endpoint primes present in row must be a
    # multiple of scale. Enumerate its positive partitions without CRT.
    required = [p for p,_ in factor(n)+factor(n-1) if valuation(n,row,p)]
    scale = prod(w for p,w in factor(n) if p in required)
    target = n//scale
    for a in range(1,target//3+1):
        for b in range(a,(target-a)//2+1):
            second = tuple(scale*x for x in (a,b,target-a-b))
            if all(valuation(n,second,p) == 0 for p in required):
                raise ArithmeticError((n,row,second))
    return target


def run(path, cross):
    rows = [tuple(map(int,line.split())) for line in Path(path).read_text().splitlines()]
    assert len(rows) == len(set(rows))
    if cross:
        subset = {r for r in rows if r[2] <= 257 and r[6] <= 8}
        assert subset == direct_box(257,8)
        print("independent small-box rows",len(subset),flush=True)
    weights = set()
    stats = Counter()
    small = set()
    for d,e,f,a,b,c,v in rows:
        t = d*e*f
        u = a+b+c
        n = u*v
        powers = (d,e,f)
        assert 1 < d < e < f and min(a,b,c) > 0 and v >= 2
        assert all(len(factor(m)) == 1 for m in powers)
        bases = [factor(m)[0][0] for m in powers]
        assert len(set(bases)) == 3 and (n-1) % t == 0 and t > u and gcd(a,b,c) == 1
        row = tuple(sorted((a*v,b*v,c*v)))
        assert all(valuation(n,row,p) == 0 for p in bases)
        g, avoided = endpoint_gcd(n,row)
        weights.add(n)
        stats['gcd below weight'] += g < n
        stats['avoided product above weight'] += avoided > n
        stats['rows using a nonprime power'] += any(m != p for m,p in zip(powers,bases))
        if g < n:
            small.add((n,v,row,g,avoided,direct_companion(n,row)))
    for n in sorted(weights):
        check_full_weight(n,stats)
    print("rows",len(rows),"weights",len(weights),"above one billion",sum(n > 10**9 for n in weights),
          "largest",max(weights))
    for name,value in sorted(stats.items()):
        print(name,value)
    print("distinct rows with small endpoint gcd",len(small))
    for record in sorted(small,key=lambda x:Fraction(x[3],x[0])):
        print(*record)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("rows")
    parser.add_argument("--cross",action="store_true")
    args = parser.parse_args()
    run(args.rows,args.cross)
