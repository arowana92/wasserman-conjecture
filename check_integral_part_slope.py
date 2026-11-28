"""Check gcd bounds for a part divided by the sum of the other two."""

from fractions import Fraction
from functools import lru_cache
from math import comb, gcd

from check_direction_pair_count import divisors
from check_primitive_row_scales import powers
from check_ratio_companion_bound import profile


@lru_cache(maxsize=None)
def constant(R):
    pairs = [(d,e) for d in divisors(R) for e in divisors(R//d) if gcd(d,e)==1]
    score,d,e = min((Fraction((d+e)**3,(d*e)**2),d,e) for d,e in pairs)
    u,H = profile((d,e,R*(d+e)))
    assert Fraction(u**3,H*H)==(R+1)**3*score
    balanced = max(min(d,e) for d,e in pairs)
    assert score*score*R>=64 and score*balanced*4>=27
    return score


def audit(limit):
    comparisons = cuts = 0
    for n in range(3,limit+1):
        rows = []
        for x in range(1,n//3+1):
            for y in range(x,(n-x)//2+1):
                row = (x,y,n-x-y)
                rows.append((row,comb(n,x)*comb(n-x,y),profile(row)))
        for row,coeff,(u,H) in rows:
            for singleton in set(row):
                h = gcd(singleton,n-singleton)
                R,S = singleton//h,(n-singleton)//h
                assert Fraction(u**3 * S**3,H**2 * (R+S)**3)>=constant(R)
                cuts += 1
            for _,other,_ in rows:
                D = gcd(n,coeff,other)
                J = gcd(n-1,coeff,other)
                assert (D*n)**3 * (J*H)**2 >= 27*(n-1)**2 * u**3
                comparisons += 1
    return comparisons,cuts


if __name__=="__main__":
    print("exact gcd comparisons and slope checks",audit(55))
    for R in range(1,1025):
        score = constant(R)
        pp = powers(R)
        if len(pp)==1:
            p,_ = pp[0]
            expect = Fraction(27,4) if p==2 else Fraction(64,9) if p==3 else Fraction(8)
            assert score==expect
    print("exact optimizing witnesses for all integer slopes through 1024")
