"""Verify the two-row identities with exact coefficients, and test the power intervals."""

from collections import Counter
from itertools import combinations_with_replacement
from math import gcd, isqrt, prod

from check_endpoint_pair_bound import rows


def check_integers(limit):
    counts = Counter()
    for n in range(3, limit+1):
        m = n-1
        for (row, first), (other, second) in combinations_with_replacement(rows(n), 2):
            d = gcd(n, gcd(first, second))
            j = gcd(m, gcd(first, second))
            weight_factor = d**3*j**2
            denominator = n**3*m**2
            for a in row:
                for b in other:
                    rest = prod(row)//a * (prod(other)//b)
                    numerator = weight_factor*rest*(n*m-a*b)
                    assert numerator > 0 and numerator % denominator == 0
                    h, k = n-a, n-b
                    assert d**15*j**10*(h*k)**13 > 2**19*(n*m)**10
                    counts['rectangle choices'] += 1
                    if rest*(n*m-a*b) < denominator:
                        counts['size alone forces a prime'] += 1
                    if rest >= n*n*m and rest*(n*m-a*b) < denominator:
                        counts['new size exclusion'] += 1
            counts['row pairs'] += 1
    print(dict(counts))


def upper_offset(q, gap, degree):
    return isqrt(q**(degree+3)//(degree*gap)**3)


def check_intervals():
    counts = Counter()
    # The direct test uses the actual remainders, without replacing h by a bound.
    for degree in range(5, 18):
        for q in range(4*degree+1, 680):
            for gap in range(1, q//(4*degree)+1):
                p = q-gap
                base = q**degree
                offset = upper_offset(q, gap, degree)
                n = base + offset
                h, k = n-p**degree, offset
                assert 2*p**degree > n
                assert h*h*k*k*(h+k-1) < 16*base**2*(base-1)**2
                assert gap*(offset//(p*q**(degree//2))) <= q
                counts['upper endpoint checks'] += 1
                step = p*q**(degree//2)
                top_multiple = n//step*step
                if top_multiple >= base:
                    counts['intervals containing a multiple'] += 1
    print(dict(counts))


def fifth_power_example():
    p, q = 1_000_000_007, 1_000_000_009
    assert all(p % i and q % i for i in range(2, isqrt(q)+1))
    step = p*q*q
    upper = q**5+upper_offset(q, q-p, 5)
    n = (upper//step)*step
    h, k = n-p**5, n-q**5
    assert n % step == 0 and 0 < k < h and 2*p**5 > n
    assert h*h*k*k*(h+k-1) < 16*n**2*(n-1)**2
    assert (q-p)*(k//step) < 2*q-(n//(q*q))%q
    assert (h*k)**2 > 16*n*n*(n-1)
    print('Fifth powers:', n, 'remainders', h, k)


if __name__ == '__main__':
    check_integers(45)
    check_intervals()
    fifth_power_example()
