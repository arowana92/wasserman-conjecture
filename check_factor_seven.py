"""Exhaust the necessary monotone base-seven weights with up to 128 digits in (N-1)/7."""

from collections import defaultdict
from itertools import combinations_with_replacement
from math import comb


def repunits(digits):
    return [(7**i-1)//6 for i in range(1, digits+1)]


def three_sums(values):
    grouped = defaultdict(list)
    for a, b, c in combinations_with_replacement(values, 3):
        total = a+b+c
        grouped[total % 72].append(total)
    return grouped


def binary_survivors(digits):
    values = repunits(digits)
    odd_lengths = three_sums(values[::2])
    even_lengths = three_sums(values[1::2])
    count = 0
    survivors = []
    for residue, odd_triples in odd_lengths.items():
        even_triples = even_lengths[(11-residue) % 72]
        count += len(odd_triples)*len(even_triples)
        for a in odd_triples:
            for b in even_triples:
                t = a+b
                if not ((6*t) & (t+1)):
                    survivors.append(t)
    return count, sorted(survivors)


def small_crosscheck():
    # Generate all six repunits directly, with no split or congruence grouping.
    vals = repunits(16)
    direct = []
    count = 0
    for terms in combinations_with_replacement(vals, 6):
        t = sum(terms)
        if t % 72 == 11:
            count += 1
            if not ((6*t) & (t+1)):
                direct.append(t)
    assert (count, sorted(direct)) == binary_survivors(16)


if __name__ == "__main__":
    small_crosscheck()
    count, survivors = binary_survivors(128)
    assert count == 232664124 and survivors == [294923]
    t, = survivors
    n = 7*t+1
    # The entire trinomial is a multiple of this binomial coefficient.
    p = 3
    valuation = 0
    while p <= n:
        valuation += n//p-(6*t)//p-(t+1)//p
        p *= 3
    assert valuation == 4
    print(f"Monotone predecessors: {comb(133, 6)}. Passing modulus 72: {count}.")
    print(f"Binary survivors: {survivors}. Ternary valuation: {valuation}.")
