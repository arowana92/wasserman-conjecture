"""Find the exact minimum of sum/H over all positive directions on a plane."""

import argparse
from fractions import Fraction
from itertools import combinations_with_replacement
from math import gcd

from check_affine_endpoint import factors
from check_direction_pair_count import divisors, heavy_rows
from check_ratio_companion_bound import profile


def minimum(r,s,t):
    if min(r,s,t)<1 or gcd(r,s,t)!=1:
        raise ValueError("positive primitive coefficients required")
    best = None
    witness = None
    for d in divisors(r):
        for e in divisors(s):
            if gcd(d,e)>1:
                continue
            f = t
            while (shared:=gcd(f,d*e))>1:
                f //= shared
            q = t//f
            A,B = r//d,s//e
            # For a fixed i, choose the least positive solution for k.
            common = gcd(B,q)
            period = q//common
            for i in range(1,q+1):
                if A*i%common:
                    continue
                k = ((-(A*i//common)*pow(B//common,-1,period))%period or period)
                l = (A*i+B*k)//q
                row = (i*e*f,k*d*f,l*d*e)
                cost = Fraction(sum(row),d*e*f)
                if best is None or cost<best:
                    scale = gcd(*row)
                    best = cost
                    witness = tuple(x//scale for x in row)
    u,H = profile(witness)
    assert Fraction(u,H)==best and r*witness[0]+s*witness[1]==t*witness[2]
    return best,witness


def audit(height):
    planes = safe = no_matching = 0
    for r,s in combinations_with_replacement(range(1,height+1),2):
        for t in range(1,height+1):
            if gcd(r,s,t)>1:
                continue
            planes += 1
            value,row = minimum(r,s,t)
            safe += value>=1
            if r==s:
                f = min(Fraction(1,d)+Fraction(1,r//d)
                        for d in divisors(r) if gcd(d,r//d)==1)
                assert value==(1+Fraction(r,t))*f
            if not any(len({p,q,k})==3 for p,_ in factors(r) for q,_ in factors(s)
                       for k,_ in factors(t)):
                no_matching += 1
                assert value>1
    print("height",height,"primitive planes",planes,"minimum at least one",safe,
          "without distinct prime representatives",no_matching,flush=True)


def direct(height):
    # If u/H < 1 then u < H <= lcm(r,s,t). Enumerate the entire possible
    # heavy locus at each height, independent of the optimization routine.
    from itertools import permutations
    from math import lcm
    planes = [(r,s,t) for r,s in combinations_with_replacement(range(1,height+1),2)
              for t in range(1,height+1) if gcd(r,s,t)==1]
    cap = max(lcm(*v) for v in planes)
    rows = {row for row,_ in heavy_rows(cap)}
    cross = {v:Fraction(1) for v in planes}
    for row in rows:
        u,H = profile(row)
        for x,y,z in set(permutations(row)):
            for r,s in combinations_with_replacement(range(1,height+1),2):
                t,rem = divmod(r*x+s*y,z)
                if not rem and (r,s,t) in cross:
                    cross[r,s,t] = min(cross[r,s,t],Fraction(u,H))
    for plane,score in cross.items():
        assert score==min(1,minimum(*plane)[0])
    print("independent plane threshold checks",len(planes),"sum bound",cap,flush=True)


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plane",type=int,nargs=3,metavar=("R","S","T"))
    parser.add_argument("--height",type=int,default=40)
    args=parser.parse_args()
    if args.plane:
        print(minimum(*args.plane))
    else:
        direct(7)
        audit(args.height)
