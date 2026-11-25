"""Check normalized binomial products and the explicit inequalities in the endpoint reduction."""

from fractions import Fraction
from functools import lru_cache
from math import comb, factorial, lcm, prod

from check_endpoint_pair_bound import rows
from prime_certificates import primes_through


@lru_cache(None)
def binomial_product(x, cutoff):
    return prod(comb(x, j) for j in range(1, cutoff+1))


@lru_cache(None)
def valuation(n, p):
    exponent = 0
    while n % p == 0:
        exponent += 1
        n //= p
    return exponent


def check_prime_divisibility():
    checks = beyond_small_layers = 0
    for n in range(3, 41):
        primes = primes_through(n)
        for parts, coefficient in rows(n):
            avoided = [p for p in primes if coefficient % p]
            for h in range(1, min(n, 10)+1):
                left = binomial_product(n, h)
                small = lcm(*range(1, h+1))
                # Test each prime separately on a single row which actually avoids it.
                targets = [(p, valuation(left, p) - h*valuation(small, p))
                           for p in avoided if left % p == 0]
                for i in range(min(h, parts[0])+1):
                    for j in range(min(h-i, parts[1])+1):
                        cutoffs = (i, j, h-i-j)
                        if cutoffs[2] > parts[2]:
                            continue
                        factors = [binomial_product(x, t) for x, t in zip(parts, cutoffs)]
                        for p, needed in targets:
                            assert sum(valuation(f, p) for f in factors) >= needed
                            checks += 1
                            beyond_small_layers += needed > 0
    print(f"Exact prime and cutoff checks: {checks}, requiring the row: {beyond_small_layers}.")


def check_constants():
    h = 996
    # Taylor sum with a geometric upper bound on the omitted tail.
    e_low = sum(Fraction(1, factorial(j)) for j in range(10))
    e_high = e_low + Fraction(11, 10*factorial(10))
    assert e_high**46 < 10**20
    assert e_high**18 < 10**8 < 148306019
    assert e_low**7 > h and (2 + 7)*Fraction(1, h) < Fraction(1, 100)
    assert Fraction((h+5)**2, 2*h*h) < Fraction(51, 100)
    assert Fraction(5, 8) + Fraction(1, h) < Fraction(63, 100)
    assert Fraction(3, 4) - Fraction(1, 10**8) - Fraction(9, h) > Fraction(73, 100)
    assert (Fraction(1, 12)-Fraction(1, h))*18 > Fraction(148, 100)
    assert Fraction(104+51+63-73, 100) < Fraction(148, 100)
    print("All seven size and logarithm estimates have strict rational margins.")


if __name__ == "__main__":
    check_prime_divisibility()
    check_constants()
