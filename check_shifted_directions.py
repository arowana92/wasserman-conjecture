"""Bound the denominator of a shifted row at any scale in an endpoint counterexample."""

from argparse import ArgumentParser
from collections import Counter
from functools import lru_cache
from itertools import combinations
from math import gcd

from carry_blocks import CarryGraph, intersection
from check_carry_intersections import dimension_percent
from check_prime_power_cofactors import carry_language, shapes
from check_shifted_prefixes import can_take, floor_cover
from search_endpoint import atoms, smallest_factors


language = lru_cache(maxsize=32768)(carry_language)


def extend(residue, modulus, target, coprime_modulus):
    """Intersect with one residue at a coprime modulus."""
    digit = (target-residue)*pow(modulus, -1, coprime_modulus) % coprime_modulus
    return residue + modulus*digit, modulus*coprime_modulus


def repeated_binary_order(row):
    return len({x & -x for x in row}) < len(row)


def binomial_divisor(x, y):
    return (x+y)//gcd(x, y)


def carry_tests(m, u, slopes, other, parity, progression, first, second, flexible):
    start, step = progression
    v = (m*start+1)//u
    delta = m*(step//u)
    a, b, c = slopes
    first_divisor = binomial_divisor(a, b)
    second_divisors = [binomial_divisor(x, y) for x, y in combinations(other, 2)]

    def accepts(p, shifted):
        if shifted and first_divisor % p == 0:
            return False
        if not shifted and any(k % p == 0 for k in second_divisors):
            return False
        digits = tuple(i for i in range(p)
                       if (u % p or (v + delta*i) % p)
                       and (p != 2 or (start + step*i) % 2 == parity))
        if not digits:
            return False
        if shifted:
            return language(p, tuple(x*step for x in slopes),
                            (a*start, b*start, c*start+1), digits)
        return language(p, tuple(x*delta for x in other),
                        tuple(x*v for x in other), digits)

    return (all(accepts(p, True) for p in sorted(first))
            and all(accepts(p, False) for p in sorted(second))
            and all(accepts(p, True) or accepts(p, False) for p in sorted(flexible)))


def possible_scales(m, factors, shape, spf):
    u, a, b, c, *other = shape
    if gcd(a, b, c) != 1:
        return
    slopes = a, b, c
    initial = -pow(m, -1, u) % u, u
    residue, modulus = initial
    forced_first = {p for p, _ in atoms(u, spf)}
    forced_second = set()
    for p, w in factors:
        take_first = can_take(w, p, a, b)
        take_second = sum(x % w == 0 for x in other) >= 2
        if not take_first and not take_second:
            return
        if not take_second:
            forced_first.add(p)
            residue, modulus = extend(residue, modulus, pow(a+b, -1, w), w)
        if not take_first:
            forced_second.add(p)

    for parity in (0, 1):
        if modulus % 2 == 0 and residue % 2 != parity:
            continue
        full = extend(residue, modulus, parity, 2) if modulus % 2 else (residue, modulus)
        first, second = forced_first.copy(), forced_second.copy()
        if parity and m % 2 and u % 2:
            second.add(2)  # The full power of 2 in N lies in V.
        if parity == 0:
            if m % 2 == 0:
                two_part = m & -m
                if sum(x % (2*two_part) == 0 for x in other) < 2:
                    continue
                second.add(2)
            elif repeated_binary_order(slopes):
                second.add(2)
        if repeated_binary_order(other):
            first.add(2)
        flexible = ({2} | {p for p, _ in factors}) - first - second
        # The cheaper graph at the progression modulo u can already exclude the case.
        if all(carry_tests(m, u, slopes, other, parity, prog, first, second, flexible)
               for prog in dict.fromkeys((initial, full))):
            yield parity, full, sorted(first)


def verify(limit, digits):
    spf = smallest_factors(limit)
    counts = Counter()
    denominators = set()
    power_cases = {k: set() for k in (2, 3, 105)}
    for m in range(3, limit+1):
        factors = atoms(m, spf)
        for shape in shapes(m, factors):
            counts["preliminary shapes"] += 1
            for parity, (residue, step), bases in possible_scales(m, factors, shape, spf):
                counts["shape and parity pairs"] += 1
                denominators.add(m)
                for k, possible in power_cases.items():
                    if any(pow(v, k, step) == residue for v in range(step)):
                        possible.add(m)
                slopes = shape[1:4]
                graphs = [CarryGraph(p, slopes, (0, 0, 1)) for p in bases]
                bounds = sorted((dimension_percent(g), g.base) for g in graphs)
                if len(bounds) < 2 or bounds[0][0]+bounds[1][0] >= 100:
                    raise ArithmeticError(("two-base dimension bound inconclusive", shape, bounds))
                scales, nodes, _ = intersection(graphs, 10**digits, step, residue)
                alternate, intervals = floor_cover(slopes, bases, step, residue, 10**digits)
                if scales != alternate or scales:
                    raise ArithmeticError(("review simultaneous carry checks", shape, scales, alternate))
                counts["prefix nodes"] += nodes
                counts["factorial-floor intervals"] += intervals
                print(f"m={m} shape={shape} W={residue} mod {step} parity={parity}"
                      f" bases={bases} dimension bounds={bounds[:2]}", flush=True)
    print(dict(counts))
    print("Possible denominators:", sorted(denominators))
    for k, possible in power_cases.items():
        print(f"Possible denominators for powers with exponent {k}:", sorted(possible))


if __name__ == "__main__":
    arg = ArgumentParser(description=__doc__)
    arg.add_argument("--limit", type=int, default=1000)
    arg.add_argument("--digits", type=int, default=1000)
    args = arg.parse_args()
    if args.limit < 3 or args.digits < 1:
        arg.error("use a denominator bound of at least three and a positive digit bound")
    verify(args.limit, args.digits)
