"""Finite arithmetic checks for cubic exponent spacing and the initial window."""
from fractions import Fraction
from itertools import combinations
from math import gcd
from check_endpoint_pair_bound import prime_powers


def orders(z):
    values = dict()
    for sign, n in ((1, abs(z.numerator)), (-1, z.denominator)):
        for p, w in prime_powers(n):
            e = 0
            while w > 1:
                w //= p
                e += 1
            values[p] = sign*e
    return values


def independent_pair(x, y):
    left, right = orders(x), orders(y)
    return any(left.get(r,0)*right.get(s,0) != left.get(s,0)*right.get(r,0)
               for r,s in combinations(left.keys() | right.keys(), 2))


def test():
    directions = checks = shared = nontrivial = 0
    for B,c,p in ((5,7,37),(7,11,79),(11,17,191),(17,11,193),(13,17,223)):
        assert p>B*c and B % c not in (1,c-1) and c % B not in (1,B-1)
        for b in range(1, B//2+1):
            for d in range(1, c//2+1):
                if gcd(b,d) != 1:
                    continue
                A = b*c-B*d
                x, y = Fraction(-d,b), Fraction(-b,A)
                assert p not in orders(x) and p not in orders(y)
                assert independent_pair(x,y) == ((b,d) != (1,1))
                directions += 1
                for h in (1,2,3,5):
                    P = p**h
                    D = B*P+c
                    line = b*P+d
                    entries = []
                    for a in range(1, 45):
                        # This uses the original weight and the primitive linear form.
                        F = gcd(1+p**a*D, line)
                        assert gcd(F,p*b*d*A) == 1
                        assert (b*P+d) % F == (A*p**a+b) % F == 0
                        q = (2*a+h)//(2*h)
                        r = a-q*h
                        numerator = p**max(r,0)*(-d)**q*A+p**max(-r,0)*b**(q+1)
                        assert numerator % F == 0
                        assert numerator != 0
                        checks += 1
                        nontrivial += F>1
                        if F>1:
                            entries.append((a,F))
                    for (a,f),(ap,fp) in combinations(entries,2):
                        G = gcd(f,fp)
                        assert G*line >= f*fp and pow(p,ap-a,G) == 1 % G
                        if G>2:
                            assert (x.numerator-x.denominator*pow(p,h,G)) % G == 0
                            assert (y.numerator-y.denominator*pow(p,a,G)) % G == 0
                            shared += 1
    assert min(shared,nontrivial)>0
    assert 374**6 < 191**7
    assert 7*64*(1300000+2*10000) < 6*10000**2
    assert 7*64*(600**2//100+2*600) <= 6*600**2
    print(directions,"primitive directions,",checks,"weight congruences,",
          nontrivial,"nontrivial gcds,",shared,"pairs with shared modulus above two")


if __name__ == "__main__":
    test()
