"""Check resultant bounds, digit comparisons, and translated-cubic examples."""

from collections import Counter
from itertools import product
from math import gcd

from flint import fmpz_mat

from check_clustered_predecessors import binomial_avoids
from check_cofactor_digit_gap import atoms
from check_scaled_row_blocks import certify_prime
from check_three_factors import powers, valuation


def at(coeffs, x):
    value = 0
    for a in reversed(coeffs):
        value = value*x+a
    return value


def resultant(f, g):
    """Sylvester determinant, with coefficients supplied in increasing order."""
    g = list(g)
    while g[-1] == 0:
        g.pop()
    m, n = len(f)-1, len(g)-1
    rows = []
    for coeffs, shifts in ((f,n),(g,m)):
        for i in range(shifts):
            rows.append([0]*i + list(reversed(coeffs)) + [0]*(shifts-i-1))
    return abs(int(fmpz_mat(rows).det()))


def check_polynomials():
    counts = Counter()
    for h in range(2,8):
        for s in range(1,h):
            for d in range(1,4):
                for b,a,v in product(range(d*s*h+1),range(d*(s+h)+1),range(d+1)):
                    F = (b,a,v)
                    if not any(F):
                        continue
                    for sign in (-1,1):
                        G = (sign,d*s*h,d*(s+h),d)
                        res = resultant(G,F)
                        assert 1 <= res <= d**5*(2*h+1)**6
                        for P in (13,64,101,256):
                            assert res % gcd(at(F,P),at(G,P)) == 0
                        counts['resultant bounds'] += 1
                    # Evaluate parts at the three relevant shifted digits.
                    constant = -h*at(F,-h)
                    linear = b-2*a*h+3*v*h*h
                    quadratic = a-3*v*h
                    assert abs(constant) <= 4*d*h**3 and abs(linear) <= 8*d*h*h
                    assert abs(quadratic) <= 5*d*h
                    if constant in (0,1) and linear >= 0:
                        assert quadratic < 0
                        counts['forced negative quadratic coefficients'] += 1
    return counts


def joint_triple(n, P, p, w, r):
    total = n//(P*w)
    eligible = {i for i in range(1,total-1)
                if binomial_avoids(n, P*w*i, p) and binomial_avoids(n, P*w*i, r)}
    for x in sorted(eligible):
        if 3*x > total:
            break
        for y in range(x,(total-x)//2+1):
            z = total-x-y
            if y not in eligible or z not in eligible:
                continue
            row = tuple(P*w*i for i in (x,y,z))
            if valuation(n,row,p) == valuation(n,row,r) == 0:
                return row
    return None


def check_unshifted():
    counts = Counter()
    for p,P in atoms(180):
        for h in range(2,10):
            for s in range(1,h):
                for d in range(1,5):
                    n = d*P*(P+s)*(P+h)
                    H = n-d*P**3
                    for j in (s,h):
                        for r,w in powers(P+j).items():
                            c = (P+j)//w
                            if n//(P*w)>350:
                                continue
                            local = (H < min(w**3,P**3) and n//w**3 >= d*c**3-1
                                     and j*(H//(P*w)) <= w)
                            if not local and (P+h+s+d+j)%53:
                                continue
                            row = joint_triple(n,P,p,w,r)
                            if local:
                                assert row is None, (P,s,h,d,w,row)
                                counts['incompatible scale assignments'] += 1
                                if c>1:
                                    counts['with a proper power divisor of a translate'] += 1
                            elif row:
                                counts['compatible controls'] += 1
    return counts


def shifted_digit_examples():
    # Start with positive polynomials having a common root at the unshifted
    # node and another at the largest shift. Check the actual carries at w.
    counts = Counter()
    for p,P in atoms(1600):
        for h in range(2,9):
            for r,w in powers(P+h).items():
                c = (P+h)//w
                for s in range(1,h):
                    for d in range(3,6):
                        n = d*P*(P+s)*(P+h)+1
                        if w <= 4*max(4*d*h**3+1,8*d*c*h*h,5*d*c*c*h):
                            continue
                        # All three have small linear digits and negative
                        # quadratic coefficients in the shifted variable.
                        M = P*(P+s)*(P+h)
                        for k in range(1,d-1):
                            row = (k*M,(d-1-k)*M,M+1)
                            assert min(row)>0 and sum(row)==n and valuation(n,row,r)>0
                            assert sum(x//w**2 for x in row) == n//w**2
                            assert sum(x//w**3 for x in row) < n//w**3
                            counts['shifted rows with a forced carry'] += 1
    return counts


def large_example():
    P = 1000000000046111699137200127
    c, e = 2**20, 3**12
    parts = ((P+1)//c,(P+2)//e)
    assert c*parts[0] == P+1 and e*parts[1] == P+2
    for z in (P,*parts):
        certify_prime(z)
    assert min(c,e)**6 > P
    def sufficient(d):
        return (P**3 >= d**16*(c*e)**3*5**18 and P >= 16*d*max(c,e)**3*5
                and P >= 16*d*max(c,e)**2*25)
    assert all(sufficient(d) for d in range(1,119)) and not sufficient(119)
    print('Prime certificates and both weights at each of 118 multipliers verified.')


if __name__ == '__main__':
    print('Polynomials:', check_polynomials())
    print('Unshifted:', check_unshifted())
    print('Shifted:', shifted_digit_examples())
    large_example()
