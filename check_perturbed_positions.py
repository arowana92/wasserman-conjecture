"""Certify the seven-unit family for the 39 parameters below the general cutoff."""

from collections import Counter
from functools import lru_cache
from itertools import combinations
from math import prod

import sympy as sp


positions = tuple(range(6))
variables = sp.symbols("X1:6")
T = sp.Symbol("T")


def complement(block):
    return tuple(i for i in positions if i not in block)


def allocations():
    for x in combinations(positions, 2):
        for size in (2, 3):
            for y in combinations(complement(x), size):
                if size == 3 or x < y:
                    yield x, y


@lru_cache(maxsize=None)
def projection(first, second):
    def directional(block):
        return sum((variables[i-1] if i else 1)*T**(i-min(block)) for i in block)

    resultant = sp.resultant(directional(first), directional(second), T)
    assert resultant != 0
    terms = sp.Poly(resultant, *variables).terms()
    assert len({sum(i*c for i, c in enumerate(index, 1)) for index, _ in terms}) == 1
    return terms


@lru_cache(maxsize=None)
def bounds(first, second, d):
    # The common term in the original exponents cancels by weighted homogeneity.
    coefficients = Counter()
    for index, c in projection(first, second):
        exponent = sum(v*d**i for i, v in enumerate(index))
        coefficients[exponent] += int(c)
    terms = sorted((k, c) for k, c in coefficients.items() if c)
    assert terms, (first, second, d)
    low = terms[0][0]
    terms = [(k-low, c) for k, c in terms]
    degree = terms[-1][0]
    norm = sum(abs(c) for k, c in terms)
    if len(terms) > 1 and degree-terms[-2][0] < norm.bit_length():
        # Every positive integer root divides the nonzero constant coefficient.
        for u in sp.divisors(abs(terms[0][1])):
            if u >= 2:
                residue = sum(c*pow(int(u), k, 1009) for k, c in terms) % 1009
                assert residue or sum(c*int(u)**k for k, c in terms) != 0
    # Otherwise at every integer u >= 2 the leading term dominates the rest.
    return degree, norm


def check():
    count = 0
    for d in range(2, 41):
        span = 5*d**5+d**4
        for x, y in allocations():
            divisors = [bounds(x, y, d), bounds(x, complement(x), d),
                        bounds(y, complement(y), d)]
            powers = (3, 2, 2)
            gap = 2*span-sum(k*w for (k, c), w in zip(divisors, powers))
            norm = prod(c**w for (k, c), w in zip(divisors, powers))
            assert gap >= 0
            assert gap >= norm.bit_length() or 27*2**gap > norm
            count += 1
    assert count == 4095
    print('Seven-unit allocation and parameter pairs', count)
    print('Distinct polynomial specializations', bounds.cache_info().currsize)


if __name__ == '__main__':
    check()
