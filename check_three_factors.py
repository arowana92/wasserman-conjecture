"""Independent check of the three-factor counts at small bounds, and an explicit row."""

from itertools import combinations
from math import gcd, isqrt, prod


def powers(n):
    factors = {}
    for p in range(2, isqrt(n) + 1):
        if n % p == 0:
            factors[p] = 1
            while n % p == 0:
                n //= p
                factors[p] *= p
        if n == 1:
            break
    if n > 1:
        factors[n] = n
    return factors


def valuation(n, row, p):
    result, w = 0, p
    while w <= n:
        result += n // w - sum(a // w for a in row)
        w *= p
    return result


def check_small(bound):
    counts = [0]*5
    atoms = []
    for w in range(2, bound + 1):
        decomposition = powers(w)
        if len(decomposition) == 1:
            atoms.append(next(iter(decomposition.items())))
    for triple in combinations(atoms, 3):
        factors = dict(triple)
        if len(factors) != 3:
            continue
        m = prod(factors.values())
        counts[0] += 1
        n = m + 1
        row = [(m//w)*pow(m//w, -1, w) for w in factors.values()]
        if sum(row) != n:
            continue
        counts[1] += 1
        if gcd(*row) < 3:
            continue
        counts[2] += 1
        if any(valuation(n, row, p) for p in factors):
            continue
        counts[3] += 1
        u = prod(w for p, w in powers(n).items() if valuation(n, row, p) == 0)
        counts[4] += u >= 3
    print("Small bound", bound, "stage counts", *counts)
    return counts


def check_example():
    n = 222788354
    row = (86617810, 45491401, 90679143)
    if sum(row) != n:
        raise ArithmeticError(row)
    primes = {**powers(n), **powers(n-1)}
    values = {p: valuation(n, row, p) for p in primes}
    expected = {2: 29, 167: 4, 283: 0, 2357: 2, 29: 0, 127: 0, 241: 0, 251: 5}
    if values != expected or prod(primes.values()) != n*(n-1):
        raise ArithmeticError(values)
    avoided = prod(primes[p] for p in values if values[p] == 0)
    if avoided != 251191649 or avoided <= n:
        raise ArithmeticError(avoided)
    print("Example weight", n, "avoided product", avoided)


if __name__ == "__main__":
    if check_small(500) != [235224, 117559, 8821, 2, 0]:
        raise ArithmeticError("The independent stage counts disagree.")
    check_example()
