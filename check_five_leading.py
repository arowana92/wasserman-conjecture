"""Leading-digit certificates at cofactor five and exact projected intervals."""

import argparse
from fractions import Fraction
from functools import lru_cache
from itertools import product
from pathlib import Path


@lru_cache(maxsize=None)
def digits(p):
    return [(a,b,c) for a,b,c in product(range(p),repeat=3) if a+b+c < p]


def leading_depth(p, slopes, denominator, limit=90, max_width=500):
    # Strict caps follow from v_i < (sum(v)+1)*part_i/weight.
    front = [(0,0,0)]
    for depth in range(1,limit+1):
        following = []
        for a,b,c in front:
            for i,j,k in digits(p):
                v,w,z = p*a+i,p*b+j,p*c+k
                n = v+w+z
                if (n and v*denominator < (n+1)*slopes[0]
                      and w*denominator < (n+1)*slopes[1]
                      and z*denominator < (n+1)*slopes[2]):
                    following.append((v,w,z))
        if not following:
            return depth
        if len(following) > max_width:
            return None
        front = following
    return None


def scalar31():
    P, c = 1 << 31, 4_192_014_075
    a,b,f = 3*P,c,2*P-c
    scale = 1 << 16
    assert leading_depth(3,(a*scale,b*scale+1,f*scale+1),5*P*scale) == 28
    worst = []
    for offset in (0,1):
        max_first = 0
        for w in range(1,scale):
            parts = (a*w,b*w,f*w+offset)
            modulus, step = 3,1
            while sum(t % modulus for t in parts) < modulus:
                modulus *= 3
                step += 1
                assert modulus <= 3*sum(parts)
            max_first = max(max_first,step)
        worst.append(max_first)
    print("Large scalars excluded at 28 leading trits. Small-scalar maximum first layers:",worst)


def extended(rows):
    expected = {36:0,37:0,38:0,39:4,40:1,41:0,42:4,43:6,44:8,45:9}
    assert len(set(rows)) == len(rows) == sum(expected.values())
    assert all(sum(k==j for k,_,_ in rows)==count for j,count in expected.items())
    left = []
    proofs = []
    for k,c,q in rows:
        P = 1 << k
        assert 0 < c < 2*P and c%15 == 0 and c*q%(2*P) == 1 and 0 < q < 2*P
        for p in (2,3,5):
            depth = leading_depth(p,(3*P*q,c*q,(2*P-c)*q+1),5*P*q)
            if depth:
                proofs.append((k,c,p,depth))
                break
        else:
            left.append((k,c,q))
    assert sorted(left) == [(44,8136949755,28147172193075),(44,32985432719355,26828080362291)]
    print("Certified by (k, coefficient, base, depth):",proofs)
    print("After the separate 58-trit and short-input check the remaining new class is:",
          (44,32985432719355,26828080362291))



def narrow_prefix():
    # Remove the first of the two surviving templates at 44. The stronger
    # rational cap applies above the cutoff, and the short binary inputs can
    # be read directly below it.
    k,c,q = 44,8136949755,28147172193075
    P = 1 << k
    limit = 3**58
    slopes = (3*P*limit,c*limit,(2*P-c)*limit+1)
    assert leading_depth(3,slopes,5*P*limit,max_width=5000) == 58
    front = [q]
    for j in range(k+1,limit.bit_length()):
        bit = 1 << j
        following = []
        for x in front:
            for w in (x,x+bit):
                a,b,z = 3*P*w,c*w,(2*P-c)*w+1
                if w<limit and not ((3*w & 2*w | a & b | (a | b) & z) & bit):
                    following.append(w)
        front = following
    assert len(front) == 7
    for w in front:
        a,b,z = 3*P*w,c*w,(2*P-c)*w+1
        assert a & b or (a | b) & z
    print("At 44 the 58-trit exclusion and seven short binary prefixes remove c =",c)


GRID = 1 << 52
MAX_DENOMINATOR = 5*(1 << 20)


def projected(p, height):
    # Store (total, 5*first - 3*total, second) at each prefix.
    steps = [(a+b+c,2*a-3*(b+c),b) for a,b,c in digits(p)]
    front = [(0,0,0)]
    for _ in range(height):
        front = [(p*n+s,p*d+v,p*b+w) for n,d,b in front for s,v,w in steps
                 if -2 <= p*d+v <= 2 and p*n+s > 0]
    intervals = sorted((GRID*b//(n+1),GRID*(5*b+2+d)//(5*(n+1))) for n,d,b in front)
    merged = []
    for x,y in intervals:
        if merged and merged[-1][1] >= x:
            merged[-1][1] = max(merged[-1][1],y)
        else:
            merged.append([x,y])
    return merged


def intersection(left, right):
    i = j = 0
    result = []
    while i<len(left) and j<len(right):
        a,b = left[i]
        x,y = right[j]
        if max(a,x) <= min(b,y):
            result.append((max(a,x), min(b,y)))
        if b<y:
            i += 1
        else:
            j += 1
    return result


def interval_check():
    assert GRID > 2*MAX_DENOMINATOR**2
    common = [(0,2*GRID//5)]
    for p,h in [(2,20),(3,12),(5,8)]:
        common = intersection(common,projected(p,h))
    positive = [(x,y) for x,y in common if x<y]
    bound = sum(y-x+1 for x,y in positive)
    assert (len(common),len(positive)) == (12978,11742)
    assert bound*1_000_000 < 15382*GRID
    end = max(y for x,y in positive if x < GRID//5)
    exact = Fraction(end,GRID).limit_denominator(MAX_DENOMINATOR)
    assert exact == Fraction(95574,2286935) and exact < Fraction(209,5000)
    print("Joint cover:",len(positive),"intervals of positive length; upper length",bound,"/",GRID)
    print("Largest permitted smaller fraction:",exact)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("rows",type=Path,nargs="*",
                        help="residual rows printed by the independent low-digit verifier")
    args = parser.parse_args()
    scalar31()
    narrow_prefix()
    interval_check()
    if args.rows:
        extended([tuple(map(int,line.split())) for f in args.rows for line in f.read_text().splitlines()])
