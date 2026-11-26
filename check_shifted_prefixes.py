"""Certify the shifted proportional rows at odd binary prefixes through 4095."""

from argparse import ArgumentParser
from bisect import bisect_right
from collections import Counter
from functools import lru_cache
from math import gcd, prod

from carry_blocks import CarryGraph, intersection
from check_carry_intersections import dimension_percent
from check_prime_power_cofactors import carry_language
from check_six_bit_rays import normalized_rows
from check_six_digit_gaps import complement_by_parts
from search_endpoint import atoms, smallest_factors


def can_take(w, p, x, y):
    return (x % w == 0 and y % p != 0) or (y % w == 0 and x % p != 0)


def fits(u, powers):
    return all(2*w+1 <= u for w in powers) and 27*prod(powers)**2 <= u**3


def possible_directions(bits, factors):
    cap = (1 << bits)-1
    divisors = [[] for _ in range(cap+1)]
    for u in range(3, cap+1):
        for i in range(u, cap+1, u):
            divisors[i].append(u)
    results = set()
    for a in range(1, cap):
        b = cap ^ a
        while b > a:
            common = gcd(a, b)
            if common >= 3:
                available = cap ^ (a | b)
                c = available
                while c:
                    if (a+b+c) % 2:
                        g = gcd(common, c)
                        x, y, z = a//g, b//g, c//g
                        s = x+y+z
                        required = [w for p, w in factors[s] if not can_take(w, p, x, y)]
                        for u in divisors[common//g]:
                            if fits(u, required):
                                results.add((s, x, y, z, u))
                    c = (c-1) & available
            b = (b-1) & (cap ^ a)
    return sorted(results)


def small_direction_check(factors, bits=8):
    # Reverse the enumeration. Start with all primitive, positive integer slopes.
    cap = (1 << bits)-1
    expected = set()
    for s in range(3, cap+1, 2):
        for x in range(1, s//2):
            for y in range(x+1, s-x):
                z = s-x-y
                if gcd(gcd(x, y), z) != 1:
                    continue
                if not any(((g*x & g*y) | (g*x & g*z) | (g*y & g*z)) == 0
                           for g in range(1, cap//s+1, 2)):
                    continue
                for u in range(3, gcd(x, y)+1):
                    if x % u or y % u:
                        continue
                    required = [w for p, w in factors[s]
                                if not (max(gcd(x, w), gcd(y, w)) == w
                                        and min(gcd(x, p), gcd(y, p)) == 1)]
                    if fits(u, required):
                        expected.add((s, x, y, z, u))
    assert expected == set(possible_directions(bits, factors))
    return len(expected)


def stream_accepts(p, slopes, offsets, start, step, s, u):
    # Input t=start+step*n. Unitarity requires p not dividing (s*t+1)/u at p|u.
    v = (s*start+1)//u
    dv = s*(step//u)
    digits = [i for i in range(p) if u % p or (v+dv*i) % p]
    return carry_language(p, tuple(a*step for a in slopes),
                          tuple(a*start+b for a, b in zip(slopes, offsets)), digits)


def candidate_rows(direction, factors, counts):
    s, x, y, z, u = direction
    slopes, offsets = (x, y, z), (0, 0, 1)
    start = -pow(s, -1, u) % u
    bases = sorted({2} | {p for p, _ in factors[u]})

    @lru_cache(None)
    def first(p):
        return stream_accepts(p, slopes, offsets, start, u, s, u)

    if not all(first(p) for p in bases):
        return
    counts["after first-row base tests"] += 1
    required = [(p, w) for p, w in factors[s] if not can_take(w, p, x, y) or not first(p)]
    powers = tuple(w for _, w in required)
    if not fits(u, powers):
        return
    counts["after complementary size tests"] += 1
    complementary = {tuple(sorted(r)) for r in normalized_rows(u, powers)}
    if complementary != complement_by_parts(u, powers):
        raise ArithmeticError(("complementary row mismatch", u, powers))
    for row in sorted(complementary):
        counts["complementary assignments"] += 1
        # Every predecessor base whose two zeros cannot fit in C is forced on A.
        forced = [(p, w) for p, w in factors[s] if sum(f % w == 0 for f in row) < 2]
        residue, stride = start, u
        for p, w in forced:
            assert can_take(w, p, x, y)
            residue += stride*((pow(x+y, -1, w)-residue)*pow(stride, -1, w) % w)
            stride *= w
        # Test the small required C graphs before imposing the combined progression.
        v0 = (s*start+1)//u
        if any(not carry_language(p, tuple(f*s for f in row), tuple(f*v0 for f in row),
                                  range(p)) for p, _ in required):
            continue
        prime_list = sorted(set(bases) | {p for p, _ in forced})
        if not all(stream_accepts(p, slopes, offsets, residue, stride, s, u)
                   for p in prime_list):
            continue
        yield row, prime_list, residue, stride


def floor_cover(slopes, bases, modulus, residue, limit, offsets=(0, 0, 1)):
    """Cover the progression by intervals on each of which one factorial layer carries."""
    powers = []
    for p in bases:
        layers = []
        w = p
        while w <= sum(slopes)*limit+sum(offsets):
            layers.append(w)
            w *= p
        powers.append(layers)
    steps = 0
    survivors = []
    t = residue + max(0, (2-residue+modulus-1)//modulus)*modulus
    while t <= limit:
        parts = [a*t+b for a, b in zip(slopes, offsets)]
        n = sum(parts)
        skip = -1
        for layers in powers:
            i = bisect_right(layers, n)-1
            while i >= 0:
                w = layers[i]
                i -= 1
                if w <= (skip+1)*max(slopes):
                    break
                floors = [x//w for x in parts]
                if n//w > sum(floors):
                    # Keep all denominator floors fixed. The numerator floor can rise.
                    reach = min(((q+1)*w-1-x)//a for q, x, a in zip(floors, parts, slopes))
                    skip = max(skip, reach)
        if skip < 0:
            survivors.append(t)
            skip = 0
        t += (skip//modulus+1)*modulus
        steps += 1
    return survivors, steps


def verify(bits, digits):
    cap = (1 << bits)-1
    spf = smallest_factors(cap)
    factors = [atoms(n, spf) if n > 1 else [] for n in range(cap+1)]
    if bits >= 8:
        check_bits = min(bits, 9)
        sample = small_direction_check(factors, check_bits)
        print(f"Independent enumeration through {(1 << check_bits)-1} agrees on "
              f"{sample} preliminary directions.", flush=True)
    directions = possible_directions(bits, factors)
    counts = Counter({"initial directions and scales": len(directions)})
    for direction in directions:
        s, x, y, z, u = direction
        for row, bases, residue, stride in candidate_rows(direction, factors, counts):
            graphs = [CarryGraph(p, (x, y, z), (0, 0, 1)) for p in bases]
            bounds = sorted((dimension_percent(g), g.base) for g in graphs)
            assert bounds[0][0]+bounds[1][0] < 100
            found, nodes, direct = intersection(graphs, 10**digits, stride, residue)
            independent, steps = floor_cover((x, y, z), bases, stride, residue, 10**digits)
            if found != independent:
                raise ArithmeticError(("floor-interval mismatch", found, independent))
            if found:
                raise ArithmeticError(("carry intersection requires further checks", direction, found))
            counts["surviving assignments"] += 1
            print(f"s={s} slopes={(x,y,z)} U={u} second row={row} "
                  f"t={residue} mod {stride} bases={bases} best dimension bounds={bounds[:2]} "
                  f"intervals={nodes} direct={direct} floor intervals={steps}", flush=True)
    print(dict(counts))


if __name__ == "__main__":
    cli = ArgumentParser(description=__doc__)
    cli.add_argument("--bits", type=int, default=12)
    cli.add_argument("--digits", type=int, default=1000)
    args = cli.parse_args()
    if args.bits < 2 or args.digits < 1:
        cli.error("use at least two prefix bits and a positive decimal bound")
    verify(args.bits, args.digits)
