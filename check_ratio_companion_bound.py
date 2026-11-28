"""Exact factorial checks for the ratio/companion bound and plane estimates."""

from itertools import combinations, combinations_with_replacement
from math import comb, gcd, lcm

from check_part_ratio_divisors import sharp_constant


def profile(row):
    scale = gcd(*row)
    primitive = tuple(i//scale for i in row)
    H = lcm(*(gcd(i,j) for i,j in combinations(primitive,2)))
    return sum(primitive), H


def exact_endpoint(limit):
    checked = 0
    for n in range(3,limit+1):
        rows = [(x,y,n-x-y) for x in range(1,n//3+1)
                for y in range(x,(n-x)//2+1)]
        data = [(row,comb(n,row[0])*comb(n-row[0],row[1]),profile(row)) for row in rows]
        for first,coeff_a,_ in data:
            for i,j in set(combinations(first,2)):
                g = gcd(i,j)
                a,b = i//g,j//g
                cap = sharp_constant(a,b)
                for _,coeff_b,(u,H) in data:
                    G = gcd(n*(n-1),coeff_a,coeff_b)
                    removal = gcd(n-1,lcm(cap,H))
                    assert (G*g*removal) % (u*(n-1)) == 0
                    assert G*removal >= u*(a+b)
                    checked += 1
    return checked


def plane_checks(height,lim):
    checks = sharp = 0
    for r,s in combinations_with_replacement(range(1,height+1),2):
        for t in range(1,height+1):
            if gcd(r,s,t)!=1:
                continue
            numerator = s*t+r*t+2*r*s
            for x in range(1,lim):
                for y in range(1,lim-x):
                    z,rem = divmod(r*x+s*y,t)
                    if rem or x+y+z>lim:
                        continue
                    u,H = profile((x,y,z))
                    assert u*r*s*t >= H*numerator
                    checks += 1
            if lcm(r,s,t)==r*s*t:
                u,H = profile((s*t,r*t,2*r*s))
                assert u*r*s*t==H*numerator
                sharp += 1
    return checks, sharp


def tail_counts(cutoffs):
    counts = [0]*len(cutoffs)
    X = max(cutoffs)
    for x in range(1,X//3+1):
        for y in range(x,(X-x)//2+1):
            h = gcd(x,y)
            for z in range(y,X-x-y+1):
                if gcd(h,z)>1:
                    continue
                H = h*gcd(x,z)*gcd(y,z)
                if H > x+y+z:
                    for i,cutoff in enumerate(cutoffs):
                        counts[i] += x+y+z<=cutoff
    return dict(zip(cutoffs,counts))


if __name__ == "__main__":
    print("exact ordered row-pair/ratio checks",exact_endpoint(46))
    print("plane inequalities and sharp controls",plane_checks(16,48))
    print("unordered primitive directions with H above their sum",tail_counts((32,64,128,256)))
