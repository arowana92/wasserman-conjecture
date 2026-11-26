"""Certify two independent power blocks and the four smaller integral exponent ratios."""

from collections import Counter
from itertools import combinations
from math import gcd, lcm

from flint import fmpz_poly

from check_close_powers import prime_base
from check_polynomial_powers import assignment_bounds, coeffs, common_multiple, divisor_bound, row_shapes


X = fmpz_poly([0,1])
UPPER = X**2+2*X+2
LOWER = 2*X*(X+1)


def positive_after(poly, start):
    shifted = poly(X+start)
    return shifted[0] > 0 and all(c >= 0 for c in shifted.coeffs())


def numerical_test(row):
    n = sum(row)
    u = gcd(*row)
    s = gcd(n-1,lcm(*(gcd(a,b) for a,b in combinations(row,2))))
    if u**3*s**2 >= 27*(n-1)**2:
        raise ArithmeticError(('size bound needs a further test',n,row))


def two_blocks():
    # The zero separates the coefficient allocations of 2+2X and 2+2X+X^2.
    for shape in row_shapes((2,2,0,2,2,1)):
        lower = [X*fmpz_poly(list(f[:2])) for f in shape]
        upper = [fmpz_poly(list(f[3:])) for f in shape]
        assert sum(lower) == LOWER and sum(upper) == UPPER
        yield lower, upper


def two_block_bounds(lower, upper):
    bounds = []
    for slot in range(4):
        shifted = [f-int(slot==i+1) for i,f in enumerate(lower)]
        polynomials = [upper[i]*shifted[k]-upper[k]*shifted[i] for i,k in combinations(range(3),2)]
        polynomials += [b for a,b in zip(shifted,upper) if not a]
        polynomials += [a for a,b in zip(shifted,upper) if not b]
        reduced = []
        for f in polynomials:
            if not f:
                continue
            while not f[0]:
                f = f//X
            reduced.append(coeffs(f))
        assert reduced
        bounds.append(divisor_bound(*reduced))
    # The full base-prime power from the weight is at most 2P.
    u = 2*X*fmpz_poly(list(bounds[0]))
    s = fmpz_poly(list(common_multiple(tuple(bounds[1:]))))
    return u, s


def independent_exponents():
    counts = Counter()
    minimum_m = UPPER*X**8+LOWER-1
    for lower,upper in two_blocks():
        u,s = two_block_bounds(lower,upper)
        for power in (2,3,4,5,8,25):
            for block in (power**8, prime_base(power)*power**8, power**10):
                row = [int(h(power))*block+int(l(power)) for l,h in zip(lower,upper)]
                n = sum(row)
                pair_lcm = lcm(*(gcd(x,y) for x,y in combinations(row,2)))
                assert int(u(power)) % gcd(*row) == 0
                assert int(s(power)) % gcd(n-1,pair_lcm) == 0
                counts['direct divisor checks'] += 1
        diff = 729*minimum_m**4-u**6*s**4
        assert positive_after(diff,8)
        counts['two-block shapes'] += 1
        for power in range(2,8):
            p = prime_base(power)
            if not p:
                continue
            a,b = abs(int(u(power))),abs(int(s(power)))
            if not a or not b:
                counts['exact minor specializations'] += 1
                continue  # separately excluded in the proof, including at larger P
            block = power**8
            while 27*(int(UPPER(power))*block+int(LOWER(power))-1)**2 <= a**3*b**2:
                row = [int(h(power))*block+int(l(power)) for l,h in zip(lower,upper)]
                numerical_test(row)
                counts['small independent specializations'] += 1
                block *= p
    return counts


def integral_exponents():
    counts = Counter()
    for k in range(4,8):
        g = coeffs((UPPER*X**k+LOWER)//X)
        m = X*fmpz_poly(list(g))-1
        for shape in row_shapes(g):
            a,b = assignment_bounds(g,shape,None)
            u = X*fmpz_poly(list(a))
            s = fmpz_poly(list(b))
            assert positive_after(729*m**4-u**6*s**4,34)
            counts['single-variable shapes'] += 1
            # Superset including all integer bases, regardless of primality.
            for power in range(2,34):
                row = [power*int(fmpz_poly(list(f))(power)) for f in shape]
                numerical_test(row)
                counts['small integral specializations'] += 1
    return counts


def elementary_checks():
    # The determinant specialization and quadratic discriminants in the proof.
    assert 2*UPPER-(LOWER-1) == 2*X+5
    assert 4*UPPER-(2*X+5)*(2*X-1) == 13
    assert (X+1)**2-(X+2)*2*X == 2-(X+1)**2
    assert (X+1)**2-(X+2)*(2*X-1) == 3-X-X**2

    # Enumerate all pairings of four out of the five original binary units.
    checks = 0
    for a in range(1,37):
        for j in range(1,73):
            P, Z = 2**a, 2**j
            n = (P**2+2*P+2)*Z+2*P*(P+1)
            exponents = [a+1,2*a+1,j+1,j+a+1,j+2*a]
            assert n == sum(2**i for i in exponents)
            u = n & -n
            if a == j:
                assert n == P*(P+2)**2 and bin(n).count('1') <= 3
                assert n == 32 or 27*(n-1)**2 > u**3
                continue
            assert u == 2*min(P,Z)
            for four in combinations(exponents,4):
                for k in range(1,4):
                    tail = [four[i] for i in range(1,4) if i != k]
                    gap = min(abs(four[0]-four[k]), abs(tail[0]-tail[1]))
                    assert gap <= (j if j > a else 2*a-j)
                    assert 27*(n-1)**2 > u**3*(1+2**gap)**2
                    checks += 1
    return checks


if __name__ == '__main__':
    print(independent_exponents())
    print(integral_exponents())
    print('Binary pairings:',elementary_checks())
