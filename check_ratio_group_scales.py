"""Check the two height-twelve families using digit graphs and factorial floors."""

import argparse
from bisect import bisect_right

from carry_blocks import CarryGraph, intersection
from check_affine_endpoint import factors


# Each line gives the normalized companion, the forced bases there, and
# the congruence on its multiplier V from the forced primes of the ratio row.
CASES = (((63,70,540),(3,5,7),609840,399757), ((11,44,112),(2,11),174636,141455))


def floor_scan(row,bases,bound,modulus,residue):
    u = sum(row)
    powers = []
    for p in bases:
        levels = [p]
        while levels[-1]<bound:
            levels.append(levels[-1]*p)
        powers.append(levels)
    v = residue
    intervals = 0
    survivors = []
    while u*v<=bound:
        n = u*v
        parts = [x*v for x in row]
        lengths = [bisect_right(levels,n) for levels in powers]
        stop = v
        for depth in range(1,max(lengths,default=0)+1):
            for levels,length in zip(powers,lengths):
                if depth>length:
                    continue
                q = levels[length-depth]
                quotients = [x//q for x in parts]
                if n//q>sum(quotients):
                    # The carry persists until one denominator floor increases.
                    change = min(((k+1)*q+x-1)//x for k,x in zip(quotients,row))
                    stop = max(stop,change)
            if stop>v:
                break
        if stop==v:
            survivors.append(v)
            stop = v+1
        v = stop+(residue-stop)%modulus
        intervals += 1
    return intervals,survivors


def run(digits):
    bound = 10**digits
    for (a,b,d),(row,bases,modulus,residue) in zip(((4,11,13860),(7,9,2772)),CASES):
        u = sum(row)
        m,rem = divmod(u*residue-1,d)
        assert rem==0
        product = d
        for p,w in factors(d):
            on_A = a%w==0 or b%w==0
            on_B = sum(x%w==0 for x in row)>=2
            assert on_A!=on_B and on_B==(p in bases)
            if on_A:
                product *= w
                assert ((a if a%p else b)*u*m)%w==1
        assert product==modulus and 0<residue<modulus
        graphs = [CarryGraph(p,row,(0,0,0)) for p in bases]
        # A block dimension bound on the first two required bases suffices.
        for graph in graphs[:2]:
            assert graph.block_bound(256)**2 < graph.base**256
        points,nodes,_ = intersection(graphs,bound//sum(row),modulus,residue)
        intervals,check = floor_scan(row,bases,bound,modulus,residue)
        assert points==check
        print(row,"decimal exponent",digits,"prefix intervals",nodes,"factorial-floor intervals",
              intervals,"joint multipliers",len(check),flush=True)
        if digits<=10000:
            assert not points


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--digits",type=int,default=10000)
    run(parser.parse_args().digits)
