"""Recover weights by CRT from a second, factorial-floor scale enumeration."""

import argparse
import subprocess
from itertools import combinations, combinations_with_replacement, permutations
from pathlib import Path
from math import gcd, lcm, prod

from flint import fmpz

from check_line_divisors import valuation


# The left ratios are exactly those surviving the predecessor divisor bound.
PAIRS = {(2,5):[(1,3)], (3,5):[(1,2),(1,4)],
         (4,5):[(1,1),(1,3),(2,3)], (5,6):[(1,1),(1,4)]}


def factor(n):
    ff=[(int(p),int(e)) for p,e in fmpz(n).factor()]
    assert all(fmpz(p).is_prime() for p,e in ff) and prod(p**e for p,e in ff)==n
    return ff


def roots(h):
    values=[0]
    for p,e in factor(h):
        q=p**e
        v=(h//q)*pow(h//q,-1,q)
        values += [(x+v)%h for x in values]
    return values


def finite_shapes():
    ratios=[pair for pair in combinations_with_replacement(range(1,7),2) if gcd(*pair)==1]
    def bound(pair):
        return 6 if pair==(3,4) else prod(pair)
    possible={frozenset((a,b)) for a,b in combinations_with_replacement(ratios,2)
              if lcm(bound(a),bound(b))>sum(a)*sum(b)}
    assert possible=={frozenset((a,b)) for b,lefts in PAIRS.items() for a in lefts}
    mixed=set()
    for u in range(3,30):
        for x in range(1,u//3+1):
            for y in range(x,(u-x)//2+1):
                z=u-x-y
                for a,b in ratios:
                    cuts=[(v//gcd(v,w),w//gcd(v,w)) for v,w in combinations((x,y,z),2)]
                    if any(lcm(a,b,v,w)<=(a+b)*(v+w) for v,w in cuts):
                        continue
                    for i,j,k in permutations((x,y,z)):
                        for r,s in combinations_with_replacement(range(1,7),2):
                            t,rem=divmod(r*i+s*j,k)
                            if rem or not 1<=t<=6 or gcd(r,s,t)!=1:
                                continue
                            for d in range((a+b)*u+1,lcm(a,b,r,s,t)+1):
                                if lcm(a,b,r,s,t)%d==0 and gcd(d,u)==1:
                                    mixed.add((a,b,d,u,x,y,z))
    assert mixed=={(3,5,60,7,1,2,4)}


def check(exe, limit):
    finite_shapes()
    for right,lefts in PAIRS.items():
        high=(limit-1)//sum(right)
        first=3 if right[0] in (3,5) else 2
        cmd=[str(Path(exe).resolve()),*map(str,right),str(high),str(first)]
        scales=pairs=extra=0
        with subprocess.Popen(cmd,stdout=subprocess.PIPE,text=True) as proc:
            for line in proc.stdout:
                h=int(line)
                scales+=1
                residues=roots(h)
                for left in lefts:
                    cap=lcm(*left,*right)
                    if gcd(h,cap)!=1:
                        continue
                    for root in residues:
                        n=1+cap*((root-1)*pow(cap,-1,h)%h)
                        if not sum(right)*h<n<=limit:
                            continue
                        g,rem=divmod(n*(n-1),cap*h)
                        if rem or not 0<sum(left)*g<n:
                            continue
                        pairs+=1
                        a=(left[0]*g,left[1]*g,n-sum(left)*g)
                        b=(right[0]*h,right[1]*h,n-sum(right)*h)
                        if any(n*(n-1)%p==0 and valuation(a,p) and valuation(b,p)
                               for p in (2,3,5)):
                            continue
                        extra+=1
                        pp={p for v in (n,n-1) for p,e in factor(v)}
                        if not any(valuation(a,p) and valuation(b,p) for p in pp):
                            raise ArithmeticError((n,a,b))
        assert proc.returncode==0
        print(right,f"scales={scales} row_pairs={pairs} further_factorizations={extra}",flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable")
    parser.add_argument("--limit",type=int,default=10**18)
    args=parser.parse_args()
    check(args.executable,args.limit)
