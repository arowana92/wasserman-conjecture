"""Close the initial intervals for cubic weights built from P, P+1 and mP+h."""

from argparse import ArgumentParser
from collections import Counter
from functools import lru_cache
from math import prod

from search_endpoint import find_row, splits
from search_even_endpoint import prime_powers
from check_close_powers import INITIAL_LIMIT


@lru_cache(None)
def factors(n):
    return tuple(prime_powers(n).items())


def inspect(n, predecessor, counts):
    if len(predecessor) < 3 or n <= INITIAL_LIMIT:
        return
    counts["weights checked beyond initial interval"] += 1
    numerator = factors(n)
    for u,first,second in splits(numerator):
        v = n//u
        if min(u,v)<3:
            continue
        for s,selected,other in splits(predecessor):
            if len(selected)<3 or s<v:
                continue
            t = (n-1)//s
            if u<2*max((w for _,w in other),default=1)+1 or 27*t*t>u**3:
                continue
            if v<2*max(w for _,w in selected)+1 or 27*s*s>v**3:
                continue
            counts["assignments"] += 1
            row = find_row(n,u,first,s,selected)
            if row is not None:
                counts["first rows"] += 1
                partner = find_row(n,v,second,t,other)
                if partner is not None:
                    raise ArithmeticError(("check potential counterexample",n,row,partner))


def verify(side):
    counts = Counter()
    checked_weights = set()
    largest = 0
    for d in (1,2,3,4,6,8,12,16):
        for m in range(1,side+1):
            for h in range(1,side+1):
                counts["families"] += 1
                B = 1+2*d*(m+h)
                limit = B**5
                for e in range(1, (limit+1).bit_length()+1):
                    for P in (2**e,2**e-1):
                        if not 2<=P<=limit:
                            continue
                        pair = factors(P),factors(P+1)
                        if any(len(x)!=1 for x in pair):
                            continue
                        third = factors(m*P+h)
                        if len(third)!=1:
                            continue
                        counts["prime-power specializations"] += 1
                        bases = {x[0][0] for x in (*pair,third)}
                        if len(bases)<3:
                            continue
                        n = 1+d*P*(P+1)*(m*P+h)
                        if n in checked_weights:
                            continue
                        checked_weights.add(n)
                        largest = max(largest,n)
                        pred = {}
                        for block in (factors(d),*pair,third):
                            for prime,w in block:
                                pred[prime] = pred.get(prime,1)*w
                        assert prod(pred.values())==n-1
                        inspect(n,tuple(pred.items()),counts)
    print(f"Box through {side}: {dict(counts)}, {len(checked_weights)} distinct weights, "
          f"largest {largest}.")


if __name__ == "__main__":
    cli = ArgumentParser(description=__doc__)
    cli.add_argument("side",nargs="?",type=int,default=32)
    verify(cli.parse_args().side)
