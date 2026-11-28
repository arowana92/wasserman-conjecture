"""Compute the gcd of multinomials on r*x + s*y = t*z, for positive x,y,z."""

import argparse
from heapq import heappop, heappush
from itertools import product
from math import gcd, prod


def prime_factors(n):
    factors = []
    p = 2
    while p * p <= n:
        if n % p == 0:
            factors.append(p)
            while n % p == 0:
                n //= p
        p = 3 if p == 2 else p + 2
    if n > 1:
        factors.append(n)
    return factors


def possible_primes(r, s, t):
    return sorted(set(prime_factors((r+t)//gcd(r,t))) |
                  set(prime_factors((s+t)//gcd(s,t))))


def valuation(row, p):
    total = sum(row)
    power = p
    value = 0
    while power <= total:
        value += total // power - sum(x // power for x in row)
        power *= p
    return value


def minimum(r, s, t, p):
    """Return the exact minimum valuation and one positive witness row."""
    common = gcd(r, s, t)
    normal = (r//common, s//common, -t//common)
    unit = next(i for i in range(3) if normal[i] % p)
    a, b = [i for i in range(3) if i != unit]
    inverse = pow(normal[unit], -1, p)
    digits_by_carry = {}

    def digit_choices(carry):
        if carry not in digits_by_carry:
            choices = []
            for x, y in product(range(p), repeat=2):
                partial = carry + normal[a]*x + normal[b]*y
                z = -inverse * partial % p
                digits = [0, 0, 0]
                digits[a], digits[b], digits[unit] = x, y, z
                next_carry = (partial + normal[unit]*z) // p
                mask = sum(1 << i for i, d in enumerate(digits) if d)
                choices.append((next_carry, x+y+z, mask, tuple(digits)))
            digits_by_carry[carry] = choices
        return digits_by_carry[carry]

    start = (0, 0, 0)
    goal = (0, 0, 7)
    parent = {}
    distance = {start: 0}
    queue = [(0, start)]
    while queue:
        cost, state = heappop(queue)
        if cost != distance[state]:
            continue
        if state == goal:
            row = [0, 0, 0]
            while state != start:
                state, digits = parent[state]
                row = [p*v+d for v, d in zip(row, digits)]
            assert min(row) > 0 and r*row[0]+s*row[1] == t*row[2]
            assert valuation(row, p) == cost
            return cost, tuple(row)
        carry, overflow, used = state
        for next_carry, digit_sum, mask, digits in digit_choices(carry):
            outgoing = (overflow + digit_sum)//p
            dest = (next_carry, outgoing, used | mask)
            candidate = cost + outgoing
            if dest not in distance or candidate < distance[dest]:
                distance[dest] = candidate
                parent[dest] = state, digits
                heappush(queue, (candidate, dest))
    raise ArithmeticError("positive row unreachable")


def line_gcd(r, s, t):
    if min(r, s, t) <= 0:
        raise ValueError("normal magnitudes must be positive")
    witnesses = {p: minimum(r, s, t, p) for p in possible_primes(r, s, t)}
    return prod(p**e for p, (e, row) in witnesses.items()), witnesses


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("r", "s", "t"):
        parser.add_argument(name, type=int)
    args = parser.parse_args()
    divisor, certificates = line_gcd(args.r, args.s, args.t)
    for prime, (exponent, witness) in certificates.items():
        print(f"p={prime}  minimum={exponent}  witness={witness}")
    print(f"gcd={divisor}")
