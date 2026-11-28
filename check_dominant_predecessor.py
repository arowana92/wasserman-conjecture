"""Check when avoided predecessor powers can reach the primitive sum of a row."""

import argparse
import subprocess
from bisect import bisect_right
from math import log
from pathlib import Path

from carry_blocks import CarryGraph
from check_direction_pair_count import heavy_rows, powers
from check_ratio_group_scales import floor_scan
from check_ratio_groups_cross import accepts


def local_shapes(limit):
    surviving = set()
    tested = 0
    for row,t in heavy_rows(limit):
        tested += 1
        u = sum(row)
        if all(accepts(p,row,(0,0,0),w,pow(u,-1,w)) for p,w in powers(t)):
            surviving.add((*row,t))
    return tested,surviving


def dimension_bound(row,p):
    length = 64
    count = CarryGraph(p,row,(0,0,0)).block_bound(length)
    guess = int(1000*log(count,p)/length)+1
    assert count**1000 < p**(length*guess)
    return guess


def descending(row,bases,bound,t,residue):
    u = sum(row)
    levels = []
    for p in bases:
        seq = [p]
        while seq[-1]<bound:
            seq.append(p*seq[-1])
        levels.append(seq)
    v = bound//u
    v -= (v-residue)%t
    points = []
    intervals = 0
    while v>0:
        n = u*v
        parts = tuple(x*v for x in row)
        lengths = [bisect_right(seq,n) for seq in levels]
        low = None
        for depth in range(1,max(lengths,default=0)+1):
            for seq,k in zip(levels,lengths):
                if depth>k:
                    continue
                q = seq[k-depth]
                quotient = n//q
                if sum(a//q for a in parts)<quotient:
                    # Decreasing the multiplier can only decrease denominator floors.
                    start = (quotient*q+u-1)//u
                    low = start if low is None else min(low,start)
            if low is not None:
                break
        if low is None:
            points.append(v)
            low = v
        v = low-1-(low-1-residue)%t
        intervals += 1
    return intervals,sorted(points)


def run(limit,digits,executable):
    tested,shapes = local_shapes(limit)
    print("primitive sum bound",limit,"group assignments",tested,"local survivors",len(shapes),
          "distinct directions",len({case[:3] for case in shapes}),flush=True)
    if executable:
        cmd = [str(Path(executable).resolve()),str(limit)]
        data = subprocess.check_output(cmd,text=True)
        assert {tuple(map(int,line.split())) for line in data.splitlines()}==shapes
    if not shapes:
        return
    first = min((sum(case[:3]),case[:3],case[3]) for case in shapes)
    maximum = forward = reverse = 0
    for x,y,z,t in sorted(shapes):
        row = (x,y,z)
        u = sum(row)
        bases = [p for p,w in powers(t)]
        combined = sum(dimension_bound(row,p) for p in bases[:2])
        maximum = max(maximum,combined)
        assert combined<1000
        r = pow(u,-1,t)
        bound = 10**digits
        left,hits1 = floor_scan(row,bases,bound,t,r)
        right,hits2 = descending(row,bases,bound,t,r)
        assert hits1==hits2
        if hits1:
            raise ArithmeticError((row,t,hits1))
        forward += left
        reverse += right
    print("first local survivor",first,"largest dimension sum numerator",maximum,
          "weight bound decimal exponent",digits,"forward intervals",forward,
          "backward intervals",reverse,flush=True)
    if limit==2000:
        assert tested==133318 and len(shapes)==713 and maximum==990
        assert first==(347,(35,130,182),455)


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--limit",type=int,default=2000)
    parser.add_argument("--digits",type=int,default=1000)
    parser.add_argument("--cross",help="compiled independent direction enumerator")
    args=parser.parse_args()
    run(args.limit,args.digits,args.cross)
