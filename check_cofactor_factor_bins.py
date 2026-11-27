"""Finite factor partitions, exact digit tails and margins for growing bases."""

from fractions import Fraction as F
from itertools import combinations_with_replacement, product
from math import gcd, prod

from check_endpoint_pair_bound import prime_powers


def check_factorizations():
    tested = bins_tested = 0
    for n in range(1, 2001):
        powers = [w for _, w in prime_powers(n)]
        divisors = [k for k in range(1, n+1) if n % k == 0]
        threshold = max(a for a, b, c in combinations_with_replacement(divisors, 3)
                        if gcd(a, b) == gcd(a, c) == gcd(b, c) == 1)
        residual = n//prod(sorted(powers, reverse=True)[:2])
        assert threshold <= residual
        if len(powers) <= 2:
            assert threshold == 1
        attained = False
        for assignment in product(range(4), repeat=len(powers)):
            # A power can stay in s or be assigned to one of the three positions.
            factors = [1, 1, 1, 1]
            for i, power in zip(assignment, powers):
                factors[i] *= power
            s, *bins = factors
            t = prod(bins)
            least = min(bins)
            assert least <= threshold
            attained |= least == threshold
            # Minimum possible normalized sum under the divisibility conditions.
            u = sum(t//d for d in bins)
            assert s*u*least >= n and s*u*threshold >= n
            bins_tested += 1
        assert attained
        tested += 1
    print(f"{tested:,} complete divisor comparisons and {bins_tested:,} assignments.")


def check_digit_tails():
    tested = 0
    for p in (2, 3, 5, 7, 11, 17, 37):
        counts = [1]
        for length in range(1, 26):
            # Convolve by 1+x+...+x^(p-1) via a sliding integer sum.
            previous = counts
            counts = []
            active = 0
            for total in range(len(previous)+p-1):
                if total < len(previous):
                    active += previous[total]
                if total >= p:
                    active -= previous[total-p]
                counts.append(active)
            assert sum(counts) == p**length
            tail = sum(num for value, num in enumerate(counts)
                       if 4*value > 3*(p-1)*length)
            # exp(-1/8) < 8/9, since exp(1/8) > 1+1/8.
            assert tail*9**length <= p**length*8**length
            tested += 1
    assert F(18, 10000) < F(1, 100)
    assert F(1, 8)-F(1, 150)-F(1, 100) > F(1, 10) > F(1, 12)
    print(f"{tested} exact digit-sum tails and both entropy margins passed.")


if __name__ == "__main__":
    check_factorizations()
    check_digit_tails()
