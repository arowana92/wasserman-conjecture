"""Integer checks for the unrestricted two-binomial row and endpoint divisors."""

from collections import Counter
from itertools import combinations
from math import comb, gcd

from check_three_factors import valuation


def digits(value, base):
    total = 0
    while value:
        value, remainder = divmod(value, base)
        total += remainder
    return total


def pair_divisor(x, y):
    return gcd(x, y)*gcd(x, y-1)*gcd(x-1, y)


def exact_coefficients():
    counts = Counter()
    for n in range(3, 151):
        endpoint = n*(n-1)
        for x in range(1, n//3 + 1):
            for y in range(x, (n-x)//2 + 1):
                z = n-x-y
                coefficient = comb(n, x)*comb(n-x, y)
                retained = endpoint//gcd(endpoint, coefficient)
                for a, b in combinations((x, y, z), 2):
                    upper = pair_divisor(a, b)
                    assert upper == gcd(a*b, a*(a-1), b*(b-1))
                    assert upper % retained == 0
                    counts['pair checks'] += 1
                    counts['sharp nontrivial'] += upper == retained > 1
    return counts


def odd_part(value):
    return value // (value & -value)


def binomial_rows():
    counts = Counter()
    for p in (2, 3, 5, 7, 11):
        pairs = [(i, d, p**i*(1+p**d)) for i in range(1, 5)
                 for d in range(1, 9)]
        # Each third part has units digit 1. Its other digits range freely.
        third_parts = [1+p*k for k in range(1, 401)]
        for index, (i, d, x) in enumerate(pairs):
            for j, e, y in pairs[index:]:
                for z in third_parts:
                    n = x+y+z
                    if valuation(n, (x, y, z), p):
                        continue
                    m = n-1
                    while m % p == 0:
                        m //= p
                    full = (n-1)//m
                    assert x % full == y % full == (z-1) % full == 0
                    counts['carry-free shapes'] += 1
                    counts['more than seven units'] += digits(n, p) > 7
                    upper = gcd(x, y, z)
                    s0, sd, se = gcd(m, x, y), gcd(m, x, z), gcd(m, y, z)
                    if p != 2:
                        # An odd weight scale >=3 forces equal two-valuations.
                        upper = odd_part(upper)
                        if upper >= 3:
                            assert x & -x == y & -y
                            assert valuation(n, (x, y, z), 2) > 0
                        s0, sd, se = map(odd_part, (s0, sd, se))
                    if d == e:
                        f = 1+p**d
                        assert sd == se == 1 and f % (upper*s0) == 0
                        assert full*m > (x+y) and m > 2*f
                        assert gcd(n*(n-1), pair_divisor(x, y)) == gcd(n*(n-1), full*f)
                        counts['equal gaps'] += 1
                    elif upper >= 3:
                        assert upper**3*(s0*sd*se)**2 < 27*m*m
                        counts['unequal with scale at least three'] += 1
                    if upper >= 3 and upper**3*(s0*sd*se)**2 >= 27*m*m:
                        assert d == e and (1+p**d)**3 >= 27*m*m
                        counts['passes companion size'] += 1
    return counts


def two_base_pairs():
    counts = Counter()
    for p in (2, 3, 5, 7, 11, 13):
        for q in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31):
            if p == q:
                continue
            for d in range(1, 81):
                f = p**d+1
                if f % q:
                    continue
                exponents = [i for i in range(1, 81) if digits(p**i*f, q) == 2]
                for i, j in combinations(exponents, 2):
                    assert (p, q, d, i, j) in ((2, 3, 1, 1, 2), (2, 3, 3, 1, 2))
                    counts['distinct lower exponents'] += 1
                counts['individual two-unit parts'] += len(exponents)
    return counts


if __name__ == '__main__':
    print(dict(exact_coefficients()), flush=True)
    print(dict(binomial_rows()), flush=True)
    print(dict(two_base_pairs()), flush=True)
