"""Check each surviving signed row independently by compulsory primes and CRT."""

from collections import Counter
from itertools import product
from math import prod

from flint import fmpz


def valuation(row, p):
    total = sum(row)
    power = p
    value = 0
    while power <= total:
        value += total//power - sum(x//power for x in row)
        power *= p
    return value


def compulsory(n, row):
    factors = [(int(p),int(p)**e) for p,e in fmpz(n).factor()]
    assert prod(v for _,v in factors)==n and all(fmpz(p).is_prime() for p,_ in factors)
    return [(p,v) for p,v in factors if valuation(row,p)>0]


def check_rows(rows):
    counts = Counter()
    for row in rows:
        n = sum(row)
        wt = compulsory(n,row)
        pred = compulsory(n-1,row)
        scale = prod(v for _,v in wt)
        total = n//scale
        mod = prod(v for _,v in pred)
        if any(2*v >= total for _,v in pred) or 27*mod**2 > total**3:
            counts["excluded by compulsory powers"] += 1
            continue
        # Normalized parts sum to total. At each predecessor power the exceptional
        # part equals total modulo that power and its two companions vanish.
        units = [total*(mod//v)*pow(mod//v,-1,v) for _,v in pred]
        for bins in product(range(3),repeat=len(units)):
            if bins and bins[0]!=0:
                continue
            parts = [sum(u for j,u in zip(bins,units) if j==i)%mod or mod for i in range(3)]
            extra = total-sum(parts)
            for i in range(0,extra+1,mod):
                for j in range(0,extra-i+1,mod):
                    other = (scale*(parts[0]+i),scale*(parts[1]+j),scale*(parts[2]+extra-i-j))
                    counts["candidate companions"] += 1
                    if all(valuation(other,p)==0 for p,_ in wt+pred):
                        raise ArithmeticError(("potential endpoint counterexample",row,other))
        counts["excluded by residue enumeration"] += 1
    return counts
