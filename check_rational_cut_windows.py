"""Check the rational-cut divisibilities against exact binomial gcds."""
from collections import Counter
from itertools import product
from math import comb, gcd


def mixed_form(n, r, s, u, t, a, b):
    if s*a != r*n:
        return False
    if t*b == u*(n-1):
        assert gcd(s, t) == 1 and u % s == 0
        return True
    if t*(b-1) == u*(n-1):
        assert gcd(s, t) == 1 and (t-u) % s == 0
        return True
    return False


def retained_controls(limit):
    tests = partial = 0
    for n in range(3,limit+1):
        coefficients = [comb(n,i) for i in range(n)]
        for a,b in product(range(1,n), repeat=2):
            common = gcd(coefficients[a],coefficients[b])
            D, J = gcd(n,common), gcd(n-1,common)
            for r,s,u,t in ((1,3,2,5),(3,7,2,9)):
                x, y = s*a-r*n, t*b-u*n
                assert D*x*y % n == 0
                assert J*(x+r)*(x+r-s)*(y+u)*(y+u-t) % (n-1) == 0
                tests += 1
                partial += 1 < D < n or 1 < J < n-1
    print(f"{tests} retained-power tests, including {partial} with a partial endpoint.")


def verify(limit):
    fractions = [(r,s) for s in (2,3,4,5,7,9) for r in range(1,s) if gcd(r,s) == 1]
    counts = Counter()
    for n in range(3,limit+1):
        m = n-1
        divisors = [gcd(comb(n,i),n*m) for i in range(n)]
        for a,b in product(range(1,n), repeat=2):
            if gcd(divisors[a],divisors[b]) != 1:
                continue
            counts["ordered cut pairs"] += 1
            for (r,s),(u,t) in product(fractions, repeat=2):
                x, y = s*a-r*n, t*b-u*n
                A, B = 2*r-s, 2*u-t
                R, T = r*(s-r), u*(t-u)
                k, rem = divmod(x*y,n)
                assert not rem
                F, G = (x+r)*(x+r-s), (y+u)*(y+u-t)
                assert F*G % m == 0
                rhs = -T*x*x-R*y*y+(B*k-A*T)*x+(A*k-R*B)*y+k*k+A*B*k+R*T
                h, rem = divmod(rhs,m)
                assert not rem
                if k:
                    assert (k*(T*x*x+R*y*y+(A*T-B*k)*x+(R*B-A*k)*y) + h*x*y ==
                            k*(k*k+A*B*k+R*T+h))
                    counts["nonzero-product conics"] += 1
                mixed = mixed_form(n,r,s,u,t,a,b) or mixed_form(n,u,t,r,s,b,a)
                if mixed:
                    counts["exact forms"] += 1
                elif not k:
                    counts["nonconstant quadratic cases"] += 1
                if abs(x*y) < n and R*abs(G) < m and T*abs(F) < m:
                    assert mixed
                    counts["strict window cases"] += 1
    assert all(counts[key] for key in ("ordered cut pairs", "nonzero-product conics",
                                       "exact forms", "nonconstant quadratic cases",
                                       "strict window cases"))
    assert gcd(gcd(comb(46,23),comb(46,10)),46*45) == 1
    assert mixed_form(46,1,2,2,9,23,10)
    print(dict(counts))


if __name__ == "__main__":
    retained_controls(80)
    verify(400)
