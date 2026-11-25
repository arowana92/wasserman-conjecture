"""Check CRT reciprocity by factorial floors and independently enumerate a progression box."""

from argparse import ArgumentParser
from collections import Counter
from fractions import Fraction
from itertools import combinations
from math import gcd, prod
from pathlib import Path
from subprocess import check_output

from check_cofactor_digit_gap import atoms
from check_three_factors import powers, valuation


def reciprocity_checks(bound):
    counts = Counter()
    for triple in combinations(atoms(bound),3):
        if len({p for p,_ in triple}) != 3:
            continue
        ws = [w for _,w in triple]
        m = prod(ws)
        n = m+1
        ks = [pow(m//w,-1,w) for w in ws]
        parts = [k*(m//w) for k,w in zip(ks,ws)]
        if sum(parts) != n:
            continue
        if len(set(ks)) == 1:
            assert any(valuation(n,parts,p) for p,_ in triple)
            counts["three equal"] += 1
        common = gcd(*ks)
        if common < 2:
            continue
        for r,full in powers(common).items():
            u = r
            while u <= full:
                aa = [k//u for k in ks]
                betas = [Fraction(m//w%w,w) for w in ws]
                caps = [(a*a*b)%1 for a,b in zip(aa,betas)]
                for cap,a,v,w in zip(caps,aa,parts,ws):
                    assert cap == Fraction((v//u)%u,u)+Fraction(a,u*w)
                avoids_u = valuation(n,parts,r)==0
                if avoids_u:
                    assert sum(caps)<1
                for i,j in combinations(range(3),2):
                    wi,wj = ws[i],ws[j]
                    a,b = aa[i],aa[j]
                    left = (a*a*Fraction(parts[j]//wi%wi,wi))%1
                    right = (b*b*Fraction(parts[i]//wj%wj,wj))%1
                    error = left+right-Fraction(a*b,wi*wj)
                    assert error in (0,1)
                    counts["reciprocal pairs"] += 1
                    constrained = a*a*betas[i]<1 and b*b*betas[j]<1
                    if constrained and error and sum(caps)<1:
                        assert valuation(n,parts,triple[i][0]) or valuation(n,parts,triple[j][0])
                        counts["forced carry at a predecessor base"] += 1
                    if a==b==1:
                        avoids_pair = all(valuation(n,parts,triple[x][0])==0 for x in (i,j))
                        assert not (avoids_u and avoids_pair)
                        counts["two equal the prime power"] += 1
                u *= r
    assert all(counts[k] for k in ("three equal","reciprocal pairs",
                                   "forced carry at a predecessor base","two equal the prime power")), counts
    return counts


def independent_search(limit,lo,hi,periods):
    primes = [p for p,w in atoms(limit) if p==w and p>2]
    shapes = {}
    counts = [0]*6
    for u in primes:
        if not lo<=u<=hi:
            continue
        for p,q in combinations((v for v in primes if v>u),2):
            modulus = p*q
            for R in range(q+1+(-pow(modulus,-1,u)-q-1)%u,periods*u*modulus,u):
                if gcd(R,modulus)!=1:
                    continue
                k,l = pow(q*R,-1,p),pow(p*R,-1,q)
                if k%u or l%u:
                    continue
                n = p*q*R+1
                left,middle = k*q*R,l*p*R
                row = left,middle,n-left-middle
                if row[-1]<=0:
                    continue
                assert all(v%u==0 for v in row)
                # Check by subtraction of floors at the relevant grouped digit.
                good = ((left//u)%u+(middle//u)%u < u and
                        (middle//p)%p <= (n//p)%p and (left//q)%q <= (n//q)%q)
                key = u,p,q,k,l
                assert shapes.setdefault(key,good)==good
                counts[2] += 1
                if not good or n%(u*u)==0 or valuation(n,row,u):
                    continue
                counts[3] += 1
                if not valuation(n,row,p):
                    counts[4] += 1
                    counts[5] += not valuation(n,row,q)
    counts[:2] = len(shapes),sum(shapes.values())
    return counts


if __name__ == "__main__":
    cli = ArgumentParser(description=__doc__)
    cli.add_argument("binary",type=Path)
    args = cli.parse_args()
    counts = reciprocity_checks(250)
    print("Reciprocity checks",dict(counts))
    # This row really does avoid both indicated predecessor primes.
    n = 18801054
    row = [17969148,446463,385443]
    assert n-1==113*379*439 and sum(row)==n
    assert row[1]==9*113*439 and row[2]==9*113*379
    values = [valuation(n,row,p) for p in (379,439,3)]
    assert values[0]==values[1]==0 and values[2]>0
    print("Two-predecessor-prime control",n,values)
    for N,scale,p,q,R,parts in (
        (5480060281140,79,1777,7507,410801,(1461760592718,2710465408801,1307834279621)),
        (42564932000094920,83,23911,24767,71875487,
         (15957176632045956,21682043933992712,4925711434056252)),
    ):
        assert N==p*q*R+1==sum(parts) and N%scale==0 and N%(scale*scale)!=0
        assert all(powers(prime)=={prime:prime} for prime in (p,q,R,scale))
        assert [valuation(N,parts,r) for r in (scale,p,q,R)]==[0,0,0,1]
        assert all(x%r==int(i==j) for i,x in enumerate(parts) for j,r in enumerate((p,q,R)))
    print("Two CRT rows checked with a single carry at their last predecessor prime.")
    box = 47,5,11,2
    expected = independent_search(*box)
    out = check_output([str(args.binary.resolve()),*map(str,box)],text=True)
    assert list(map(int,out.splitlines()[-1].split()))==expected
    print("Independent progression counts",expected)
    positive = check_output([str(args.binary.resolve()),"1597","101","101","1"],text=True)
    R = 618793
    parts = (299428363563,465048310013,285993129948)
    n = 1063*1597*R+1
    assert n==sum(parts) and R==7*109*811
    assert [valuation(n,parts,p) for p in (101,1063,1597)]==[0,0,0]
    assert positive.splitlines()[0]=="candidate 101 1063 1597 "+str(R)+" "+" ".join(map(str,parts))
    assert list(map(int,positive.splitlines()[-1].split()))==[622546,52599,622534,9,1,1]
    print("Positive control for all three tested prime bases",n)
