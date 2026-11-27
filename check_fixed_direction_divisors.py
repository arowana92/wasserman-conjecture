"""Find the exact fixed prime valuations of trinomial coefficients under dilation."""

from bisect import bisect_right
from fractions import Fraction
from heapq import heappop, heappush
from itertools import combinations_with_replacement
from math import comb, gcd
from random import Random

from check_primitive_row_scales import powers


def valuation_from_floors(parts, scale, p):
    entries = [scale*x for x in parts]
    n = sum(entries)
    total = 0
    while n:
        n //= p
        entries = [a//p for a in entries]
        total += n-sum(entries)
    return total


def minimum_valuation(parts, p):
    slopes = (*parts, sum(parts))
    zero = (0,)*len(slopes)
    best = {}
    parent = {}
    pending = []

    def relax(state, digits, cost, length):
        for d in digits:
            successor = tuple((s+d*x)//p for s, x in zip(state, slopes))
            next_cost = cost+successor[-1]-sum(successor[:-1])
            key = (next_cost, length+1)
            if successor not in best or key < best[successor]:
                best[successor] = key
                parent[successor] = (state if length else None, d)
                heappush(pending, (*key, successor))

    relax(zero, range(1, p), 0, 0)
    while pending:
        cost, length, state = heappop(pending)
        if best[state] != (cost, length):
            continue
        if state == zero:
            scale = 0
            while state is not None:
                state, d = parent[state]
                scale = p*scale+d
            assert length <= 2*sum(parts)-len(parts)
            assert valuation_from_floors(parts, scale, p) == cost
            return cost, scale
        relax(state, range(p), cost, length)
    raise AssertionError("Every input can be finished by appending zero digits")


def fixed_divisor(parts):
    n = sum(parts)
    divisor = 1
    for p in range(2, n+1):
        if powers(p) != [(p, p)]:
            continue
        if valuation_from_floors(parts, 1, p):
            valuation, _ = minimum_valuation(parts, p)
            divisor *= p**valuation
    return divisor


def interval_minimum(parts, p):
    """Reverse shortest distances on cells cut out by rational breakpoints."""
    n = sum(parts)
    breaks = sorted({Fraction(i, z) for z in (*parts, n) for i in range(z)})
    sample = [(a+b)/2 for a, b in zip(breaks, breaks[1:]+[1])]
    costs = [int(n*t)-sum(int(x*t) for x in parts) for t in sample]
    targets = [[bisect_right(breaks, (d+t)/p)-1 for d in range(p)] for t in sample]
    distance = [0]+[float("inf")]*(len(sample)-1)
    for _ in breaks:
        updated = [min(old, min(costs[k]+distance[k] for k in choices))
                   for old, choices in zip(distance, targets)]
        if updated == distance:
            break
        distance = updated
    return min(costs[k]+distance[k] for k in targets[0][1:])


def digit_sum(n, p):
    result = 0
    while n:
        n, digit = divmod(n, p)
        result += digit
    return result


def controls():
    comparisons = large = 0
    for parts in combinations_with_replacement(range(1, 9), 3):
        for p in (2, 3, 5, 7, 11):
            value, g = minimum_valuation(parts, p)
            assert value == interval_minimum(parts, p)
            assert all(valuation_from_floors(parts, h, p) >= value for h in range(1, 81))
            large += g > 80
            comparisons += 1
    exact = 0
    for parts in combinations_with_replacement(range(1, 6), 3):
        bound = fixed_divisor(parts)
        n = sum(parts)
        direct = 0
        for g in range(1, 31):
            direct = gcd(direct, comb(g*n, g*parts[0])*comb(g*(n-parts[0]), g*parts[1]))
            assert direct % bound == 0
            exact += 1
        for p, _ in powers(direct):
            value, g = minimum_valuation(parts, p)
            assert bound % p**value == 0 and bound % p**(value+1)
            direct = gcd(direct, comb(g*n, g*parts[0])*comb(g*(n-parts[0]), g*parts[1]))
        assert direct == bound
    assert minimum_valuation((2, 7, 12), 5)[0] == 1
    assert minimum_valuation((1, 4, 6), 2)[0] == 1
    return comparisons, large, exact


def family_controls():
    count = delayed = 0
    rng = Random(1117)
    for p in (2, 3, 5, 7, 11, 17, 31, 101):
        for a, b in ((2, 1), (5, 2), (11, 10), (p+3, p)):
            row = (a, p*b, (p-1)*(a-b))
            for g in (1, 7, 243, 2**100+1, rng.randrange(1, 10**24)):
                cost = valuation_from_floors(row, g, p)
                other = valuation_from_floors((p*a, b), g, p)
                assert (p-1)*(cost-other) == digit_sum(g*row[2], p)
                count += 1
            for j in (1, 2, 7, 19):
                g = (p**j-1)//gcd(p**j-1, row[2])
                for t in (1, 7, 12345):
                    assert valuation_from_floors(row, t*g, p) >= j
                    count += 1
        if p > 2:
            for k in (1, 2, 3, 19, 147):
                g = (p**k+1)//2
                row = tuple(g*x for x in (2, p-1, p))
                layers = [sum(row)//p**j-sum(a//p**j for a in row) for j in range(1, k+3)]
                assert layers == [0]*k+[1, 0]
                delayed += 1
    return count, delayed


def first_digits():
    checked = 0
    for p in range(2, 98):
        if powers(p) != [(p, p)]:
            continue
        for a in range(p):
            for b in range(p):
                predicted = (1+a+b) % p == 0 or max(a, b) == p-1 or a+b == p
                observed = all(d+(d*a % p)+(d*b % p) >= p for d in range(1, p))
                assert predicted == observed
                checked += 1
    return checked


if __name__ == "__main__":
    print("Graph comparisons, witnesses above 80, exact multinomials:", controls())
    print("Digit bounds and arbitrarily delayed carry controls:", family_controls())
    print("Normalized projective first-digit tests:", first_digits())
