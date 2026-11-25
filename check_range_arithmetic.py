"""Independent checks for the consecutive-weight proof, using Python integers."""

from functools import lru_cache
from itertools import combinations
from math import comb, isqrt, prod
from random import Random

from prime_certificates import primes_through
from search_endpoint import find_row

PRIMES = primes_through(31623)


def factor(n):
    factors = []
    for p in PRIMES:
        if p*p > n:
            break
        if n % p == 0:
            w = 1
            while n % p == 0:
                n //= p
                w *= p
            factors.append((p, w))
    if n > 1:
        factors.append((n, n))
    return factors


def subsets(factors):
    for length in range(len(factors)+1):
        for group in combinations(factors, length):
            yield list(group), prod(w for p, w in group)


def possible(group, product, target):
    return (target >= 3 and all(2*w < target for p, w in group)
            and 27*product**2 <= target**3)


def endpoint(first, last):
    count = 0
    for n in range(first, last+1):
        factors = factor(n-1)
        for group, shift in subsets(factors):
            if len(group) < 3:
                continue
            rest = [atom for atom in factors if atom not in group]
            for scaled, u in subsets(factor(n)):
                v = n//u
                if shift < v or not possible(group, shift, v):
                    continue
                if not possible(rest, (n-1)//shift, u):
                    continue
                count += 1
                if find_row(n, u, scaled, shift, group):
                    raise ArithmeticError((n, u, shift))
    return count


@lru_cache(None)
def weighted(x, t):
    return prod((x-j)**(t-j) for j in range(t))


def layer_tests():
    rng = Random(785)
    tested = free = 0
    for n in range(12, 121):
        for repeat in range(75):
            a, e = sorted(rng.sample(range(1, n), 2))
            row = a, e-a, n-e
            value = comb(n, a)*comb(n-a, e-a)
            cuts = tuple(rng.randrange(x+1) for x in row)
            height = sum(cuts)
            cost = prod(weighted(x, t) for x, t in zip(row, cuts))
            for p in primes_through(n):
                power, exponent = p, 0
                while power <= n:
                    residue = n % power
                    if value % p == 0:
                        # Largest feasible sum of three residues with this residue class.
                        residue = max(range(residue, 3*power-2, power))
                    exponent += sum(range(height-residue, 0, -3*power))
                    power *= p
                if cost % (p**exponent):
                    raise ArithmeticError((n, row, cuts, p, exponent))
                tested += 1
                free += value % p != 0
    return tested, free


def check_log_bounds():
    for e in range(30):
        for m in range(128, 256):
            smallest = (m*2**e+127)//128
            largest = ((m+1)*2**e-1)//128
            if smallest > largest:
                continue
            low = 8*e + (m**8).bit_length()-1-56
            if 2**low > smallest**8 or largest**8 > 2**(low+2):
                raise ArithmeticError((e, m))


def main():
    check_log_bounds()
    print("Weighted valuation checks, of which carry-free:", *layer_tests())
    for first, last in [(3, 20000), (99999001, 100000000)]:
        print(first, last, "size-feasible assignments:", endpoint(first, last))
    n, u, s = 222788354, 283, 887603
    row = find_row(n, u, factor(u), s, factor(s))
    if row != (45491401, 86617810, 90679143):
        raise ArithmeticError(row)
    print("Positive heavy-row check:", n, u, s, *row)


if __name__ == "__main__":
    main()
