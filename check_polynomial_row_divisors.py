"""Exact controls for polynomial divisors, low digit sums, and power-map counting."""

from collections import Counter
from itertools import combinations_with_replacement
from math import comb, gcd, lcm, prod

from sympy import Poly, QQ, ZZ, symbols

X = symbols("X")


def bezout_bound(polys):
    """Return the primitive common factor and a certified integer multiplier."""
    nonzero = [f for f in polys if not f.is_zero]
    common = nonzero[0]
    for f in nonzero[1:]:
        common = common.gcd(f)
    common = common.primitive()[1]
    if common.LC() < 0:
        common = -common
    quotients = [f.exquo(common).set_domain(QQ) for f in nonzero]
    coefficients = [Poly(1,X,domain=QQ)]
    accumulator = quotients[0]
    for f in quotients[1:]:
        left,right,accumulator = accumulator.gcdex(f)
        coefficients = [c*left for c in coefficients] + [right]
    coefficients = [c.mul_ground(1/accumulator.LC()) for c in coefficients]
    assert sum(a*f for a,f in zip(coefficients,quotients)) == Poly(1,X,domain=QQ)
    bound = lcm(*(int(v.q) for c in coefficients for v in c.all_coeffs()))
    return common,bound


def divisor_controls():
    scale = 1+X*(X+1)
    rows = [(a*X,b*X,c*X+offset)
            for a,b,c in combinations_with_replacement(range(1,5),3)
            for offset in range(3)]
    rows.extend([(X,1,X+1), (X+2,X**2-X+1,X**2+3),
                 (scale*X**2,scale*X*(X+1),scale*(X+1)**2),
                 ((X+1)**2,X*(X+1),X**2)])
    checked = exceptions = multiple_shifted = 0
    for row in rows:
        parts = [Poly(f,X,domain=ZZ) for f in row]
        H = sum(parts)
        factors,bounds = zip(*(bezout_bound([f-int(i==j) for i,f in enumerate(parts,1)])
                               for j in range(4)))
        if sum(f.degree()>0 for f in factors[1:]) >= 2:
            multiple_shifted += 1
        D = H.exquo(factors[0])*(H-1).exquo(prod(factors[1:]))
        for t in range(1,11):
            a,b,c = (int(f.eval(t)) for f in parts)
            n = a+b+c
            coefficient = comb(n,a)*comb(n-a,b)
            divisor = abs(int(D.eval(t)))
            assert coefficient*prod(bounds) % divisor == 0
            exceptions += coefficient % divisor != 0
            checked += 1
    assert exceptions and multiple_shifted
    print(f"{checked} exact trinomial coefficients, {exceptions} requiring a multiplier, "
          f"{multiple_shifted} rows with multiple shifted polynomial factors.")


def valuation(n,r):
    assert n
    v = 0
    while n % r == 0:
        n //= r
        v += 1
    return v


def digit_sum(n,r):
    total = 0
    while n:
        n,digit = divmod(n,r)
        total += digit
    return total


def modular_controls():
    root_tests = derivative_tests = carry_tests = kernels = 0
    for r in (2,3,5,7):
        for k in range(1,6):
            modulus = r**k
            for f in ((X-1)**3,r**2*(X**2-1)**2, r*X**2+X+1, 3*X+1):
                polynomial = Poly(f,X)
                evaluate = polynomial.eval
                degree = polynomial.degree()
                depth = (k-valuation(int(polynomial.LC()),r)+degree-1)//degree
                if depth > 0:
                    roots = {t % (r**depth) for t in range(modulus) if evaluate(t)%modulus == 0}
                    assert len(roots) <= degree
                    root_tests += 1
                # Count actual fibres inside each nonsingular class.
                derivative = polynomial.diff()
                u = min(2,k)
                fibres = Counter((t%r**u,int(evaluate(t))%modulus) for t in range(modulus)
                                 if derivative.eval(t) % r**u)
                assert max(fibres.values(),default=0) <= r**u
                derivative_tests += 1
            for a in range(1,31):
                fibres = Counter(pow(t,a,modulus) for t in range(modulus) if t%r)
                assert max(fibres.values()) <= 2*(r-1)*gcd(a,modulus)
                kernels += 1
            assert sum(gcd(a,modulus) for a in range(1,31)) <= (k+1)*30
            for t in range(1,min(modulus+1,100)):
                parts = (t*(t+1), (t-1)**2+1, t+2)
                weight = sum(parts)
                low_difference = sum(digit_sum(v%modulus,r) for v in parts)
                low_difference -= digit_sum(weight%modulus,r)
                powers = [r**j for j in range(1,k+1)]
                layers = [weight//w-sum(v//w for v in parts) for w in powers]
                assert low_difference == (r-1)*sum(layers)+layers[-1]
                carry_tests += 1
    print(f"{root_tests} root bounds, {derivative_tests} nonsingular fibre bounds, "
          f"{kernels} power-map bounds, {carry_tests} low-digit identities.")


if __name__ == "__main__":
    divisor_controls()
    modular_controls()
