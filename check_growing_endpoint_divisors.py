"""Exact checks of retained-bin lifting and carry-tolerant divisor bounds."""
from fractions import Fraction as F
from itertools import combinations_with_replacement
from math import gcd, prod

from check_binary_windows import tangent
from check_endpoint_pair_bound import allocate, prime_powers, rows


def margins():
    r, gamma, k = F(517,1000), F(1,8), F(33,50)
    h = lambda x: tangent(x,x)
    bounds = (h(r), (8*k-2+13*gamma)/5, h((r+gamma)/3)+2*(1-k)/3,
              h((r+gamma)/2)+3*gamma/4)
    assert max(bounds) < F(9993,10000)
    print("All four exponent bounds are strictly below 0.9993")


def weight_checks():
    count = partial = singleton = 0
    for n in range(3,51):
        weight, pred = tuple(prime_powers(n)), tuple(prime_powers(n-1))
        table = sorted(rows(n), key=lambda pair: pair[1] & -pair[1])
        for (a,ka),(b,kb) in combinations_with_replacement(table,2):
            full = gcd(ka,kb)
            missing_n, missing_m = gcd(n,full), gcd(n-1,full)
            u,v,p0 = allocate(weight,ka,kb)
            s,t,p1 = allocate(pred,ka,kb)
            assert missing_n*u*v == n and missing_m*s*t == n-1
            assert all(x%u==0 and x*(x-1)%s==0 for x in a)
            normalized = tuple(x//v for x in b)
            bins = sorted((t//gcd(t,x),x) for x in normalized)
            assert prod(d for d,x in bins) == t
            d,x = bins[0]
            other = t//d
            z,y = x//other,(missing_n*u-x)//d
            assert z*other+d*y == missing_n*u
            assert missing_m*s*d*d*y-z == (missing_m*s*d-v*z)*missing_n*u > 0
            residue = 1+missing_m*s*d*(pow(v*z,-1,d) if d>1 else 0)
            assert (n-residue) % (missing_m*s*d*d) == 0
            assert gcd(missing_m*s*d*d,v) == 1
            bit_sum = sum(bin(x).count("1") for x in a)-bin(n).count("1")
            assert (ka & -ka) == 1<<bit_sum
            for i,j in ((0,1),(0,2),(1,2)):
                assert (a[i]*a[j]*ka) % (n*(n-1)) == 0
            if a[0]==1:
                assert u==1 and a[1]%s==0 and a[2]%s==0
                assert 27*(n-1)**2 <= (s*missing_m)**2*missing_n**3
                singleton += 1
            partial += p0 or p1
            count += 1
    print(f"{count:,} exact pairs, including {partial:,} partial assignments and {singleton:,} with a singleton part in the first row")


def smooth_divisors():
    # For the identity under the fractional moment check integral test
    # exponents, where all coefficients and both sides are exact integers.
    for limit in (2,5,13,31):
        for e in (1,2,3):
            for n in range(1,351):
                prime_factors = [(p,w) for p,w in prime_powers(n) if p <= limit]
                target = prod(w**e for p,w in prime_factors)
                terms = {1:1}
                for p,w in prime_factors:
                    additions = {}
                    j = p
                    while j <= w:
                        for d,c in terms.items():
                            additions[d*j] = c*(j**e-(j//p)**e)
                        j *= p
                    terms.update(additions)
                assert all(n%d==0 for d in terms) and sum(terms.values()) == target
    print("4,200 exact smooth-part divisor sums")


if __name__ == "__main__":
    margins()
    weight_checks()
    smooth_divisors()
