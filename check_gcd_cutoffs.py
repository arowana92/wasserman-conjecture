"""Exact local and finite checks for gcd_cutoffs.md."""

from fractions import Fraction as F
from itertools import combinations_with_replacement
from math import comb, gcd, lcm, prod

from check_endpoint_pair_bound import rows
from check_normalized_layers import binomial_product, valuation
from prime_certificates import primes_through


def check_layer_bound():
    checks = clipped = needed = 0
    for n in range(3, 33):
        for h in range(1, min(n, 7)+1):
            left = binomial_product(n, h)
            small = lcm(*range(1, h+1))
            shifts = lcm(*range(n-h+1, n+1))
            primes = [p for p in primes_through(n) if left % p == 0]
            for (x, y, z), coefficient in rows(n):
                retained = gcd(coefficient, shifts)
                clipped += retained != coefficient
                required = {
                    p: valuation(left, p)-h*valuation(small*retained, p)
                    for p in primes
                }
                for i in range(min(h, x)+1):
                    for j in range(max(0, h-i-z), min(h-i, y)+1):
                        row_product = (binomial_product(x, i)*binomial_product(y, j)
                                       * binomial_product(z, h-i-j))
                        for p, minimum in required.items():
                            assert valuation(row_product, p) >= minimum
                            checks += 1
                            needed += minimum > 0
    assert clipped and needed
    print(f'Truncated valuation and cutoff comparisons: {checks}, requiring the row: {needed}.')


def check_weight_divisor():
    count = 0
    for n in range(3, 22):
        for triple in combinations_with_replacement(rows(n), 3):
            common = gcd(*(coefficient for _, coefficient in triple))
            assert common*prod(gcd(*parts) for parts, _ in triple) % n == 0
            count += 1
    print(f'Weight divisibility on exact triples of rows: {count}.')


def check_five_and_six():
    pairs = 0
    for n in range(7, 46):
        for (first, a), (second, b) in combinations_with_replacement(rows(n), 2):
            pairs += 1
            common = gcd(a, b)
            for h in (5, 6):
                if first[:2] == (1, 1) or second[:2] == (1, 1):
                    endpoint_gcd = gcd(n*(n-1), common)
                    assert endpoint_gcd**3 >= 27*(n-1)**2
                else:
                    shared = gcd(comb(n, h), common)
                    constant = 6**21 * (60*h)**6 * binomial_product(h-1, h-1)
                    assert constant*shared**6 >= n**(h-3)
                for parts in (first, second):
                    if parts[:2] != (1, 1):
                        cutoffs = (1, 2, 3) if parts[0] == 1 else (2, 2, 2)
                        assert all(t <= x for t, x in zip(cutoffs, parts))
                        assert prod(binomial_product(x, t) for t, x in zip(cutoffs, parts)) <= n**9
    print(f'Five- and six-shift checks on exact pairs: {pairs}.')


def check_margins():
    # In (3) H >= 60, log(N/H) >= (2/5) log N and log G <= (1/2) log N.
    assert (F(1, 12)-F(1, 120))*F(2, 5) > F(1, 120)
    assert binomial_product(4, 4) == 96 and binomial_product(5, 5) == 2500
    for x in range(1, 60):
        for y in range(x, 60):
            if x+y < 5:
                continue
            h = min(6, x+y)
            i = min(x, h)
            j = h-i
            assert 0 <= j <= y
            assert binomial_product(x, i)*binomial_product(y, j) <= 2**354
    print('The height-six cutoffs and the large-height logarithmic margin pass.')


if __name__ == '__main__':
    check_layer_bound()
    check_weight_divisor()
    check_five_and_six()
    check_margins()
