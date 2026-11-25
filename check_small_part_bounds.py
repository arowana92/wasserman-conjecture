"""Check the explicit constants in the unrestricted small-part and prime-gap proof."""

from fractions import Fraction
from functools import lru_cache
from math import prod

from prime_certificates import primes_through


def triangular(t):
    return t*(t+1)//2


def weighted(x, t):
    return prod((x-j)**(t-j) for j in range(t))


def staircase(d, step):
    if d <= 0:
        return 0
    length = (d+step-1)//step
    return length*d-step*length*(length-1)//2


@lru_cache(None)
def defect_powers(height):
    powers = []
    for p in primes_through(height):
        exponent, w = 0, p
        while w <= height:
            worst = 0
            for h in range(w):
                m = h+w*((3*w-3-h)//w)
                loss = (staircase(height-h, w)-staircase(height-h, 3*w)
                        - 2*staircase(height-m, 3*w))
                worst = max(worst, loss)
            exponent += worst
            w *= p
        powers.append((p, exponent))
    return powers


def defect(height):
    return prod(p**e for p, e in defect_powers(height))


def preliminary():
    x, s = 10**9, 2000
    if [defect(h) for h in (6, 9, 12)] != [540, 26460000, 107820266400000]:
        raise ArithmeticError("small-layer constants")
    # Normalize by the indicated power of 2 or 3 to avoid rounding.
    if 108*3**30 >= 2**25*x**5:
        raise ArithmeticError("height 12 row with minimum 3")
    tests = [
        (12, s**45*x**63, 3**60),
        (6, 4*s**6*x**19, 3**18),
        (9, 108*s**15*x**39, 3**18*2**20),
        (9, (s//2)**20*x**39, 3**18*2**20),
    ]
    for h, numerator, denominator in tests:
        if defect(h)*numerator >= denominator*(x-h+1)**triangular(h):
            raise ArithmeticError(("small sums", h))
    # After the gap reduction the lower weight is 10^20.
    x = 10**20
    for m in range(4, 1001):
        h = 6
        k = min(m, h)
        rest = h-k
        # For m = 4 or 5 put the remaining one or two cutoffs on distinct parts.
        cost = weighted(m, k)*x**rest
        if defect(h)*cost*x**18 >= 6**9*(x-5)**21:
            raise ArithmeticError(("individual part", m))


@lru_cache(None)
def lower_log(n):
    return (n**1024).bit_length()-1


def upper_log(n):
    return lower_log(n)+int((n & (n-1)) != 0)


def explicit_power():
    x, s = 10**20, 10**8
    # Zero denotes a part greater than 1000. The last column is the cutoff on N-a-b.
    cases = [(1,3,6,0), (1,0,21,5), (2,3,6,0), (2,0,9,2),
             (3,1,6,0), (3,2,6,0), (3,3,6,0), (3,0,6,0),
             (0,1,717,118), (0,2,6,0), (0,3,6,0), (0,0,9,1)]
    for a, m, h, k in cases:
        # The third row has three large parts. Each selected h is divisible by three.
        balanced = 3*triangular(h//3)
        log_upper = balanced*(upper_log(x)-lower_log(3))
        n_degree = balanced+triangular(k)
        log_upper += triangular(k)*upper_log(x)
        if m:
            u = (h-m+1)//2
            degree = triangular(u)+triangular(h-m-u)
            log_upper += upper_log(weighted(m, m))+degree*upper_log(x//2)
        else:
            degree = balanced
            log_upper += balanced*(upper_log(x)-lower_log(3))
        n_degree += degree
        if a:
            s_degree = triangular(h-k-a)
            log_upper += upper_log(weighted(a, a))+s_degree*upper_log(s)
        else:
            u = (h-k+1)//2
            s_degree = triangular(u)+triangular(h-k-u)
            log_upper += s_degree*upper_log(s//2)
        log_upper += sum(e*upper_log(p) for p, e in defect_powers(h))
        if 500*n_degree+199*s_degree > 500*triangular(h):
            raise ArithmeticError(("growth degree", a, m))
        if log_upper >= triangular(h)*lower_log(x-h+1):
            raise ArithmeticError(("constant", a, m))
    return len(cases)


def log_interval(x):
    z = (x-1)/(x+1)
    terms = 40
    low = 2*sum(z**(2*i+1)/(2*i+1) for i in range(terms))
    tail = 2*z**(2*terms+1)/((2*terms+1)*(1-z*z))
    return low, low+tail


def asymptotic_margin():
    # Rational upper bound on half the alternating residue-class sum.
    count = 256
    series = sum(Fraction(1, (3*j+1)*(3*j+2)) for j in range(count))
    defect_upper = (series+Fraction(1, 9*(count-1)))/2
    low2, _ = log_interval(Fraction(2))
    low3, _ = log_interval(Fraction(3))
    _, high_c = log_interval(Fraction(273, 100))
    margin = (61*low2-25*high_c)/144+low3/6-defect_upper
    if margin <= Fraction(1, 20000):
        raise ArithmeticError("critical constant")
    cases = [(1,3,6,0), (1,0,24,5), (2,3,6,0), (2,0,12,2),
             (3,1,6,0), (3,2,6,0), (3,3,6,0), (3,0,6,0),
             (0,2,6,0), (0,3,6,0), (0,0,12,0)]
    for a, m, h, k in cases:
        d = 3*triangular(h//3)
        other = triangular((h-m+1)//2)+triangular((h-m)//2) if m else d
        near = (triangular(h-k-a) if a else
                triangular((h-k+1)//2)+triangular((h-k)//2))
        if 5*(other+d+triangular(k))+2*near >= 5*triangular(h):
            raise ArithmeticError(("strict growth degree", a, m))


if __name__ == "__main__":
    preliminary()
    print("Excluded sums through 2000 and the 997 individual parts from 4 through 1000.")
    print("Verified all", explicit_power(), "cutoff cases for the exponent 199/500.")
    asymptotic_margin()
    print("Verified the strict degree cases and the rational margin at constant 2.73.")
