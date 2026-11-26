"""Independently check the exponent classes, all local exclusions and finite lift controls."""

from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from math import gcd, lcm

from carry_blocks import CarryGraph
from check_binary_two_factors import (Shape, finite_check, finite_progression, first_empty_level,
                                      generate, predecessor_table, terminal_possible)
from check_three_factors import powers


@lru_cache(None)
def inverse_powers(modulus):
    orbit = [1]
    while True:
        x = orbit[-1]
        x = (x + modulus*(x % 2))//2
        if x == 1:
            return {value: j for j, value in enumerate(orbit)}
        orbit.append(x)


def brute_joint(a, m, b, n):
    period = lcm(m, n)
    return next(((v, period) for v in range(a, period, m) if v % n == b), None)


def mask_shapes(bound):
    full = []
    for w in range(3, bound+1, 2):
        factors = powers(w)
        if len(factors) == 1:
            full.append((w, next(iter(factors))))
    found = set()
    for (P, p), (Q, q) in combinations(full, 2):
        if p == q:
            continue
        d = P*Q
        xs, ys = [], []
        mask = d
        while mask:
            mask = (mask-1) & d
            if mask and mask % Q == 0 and mask % P in inverse_powers(P):
                xs.append(mask)
            if mask and mask % P == 0 and mask % Q in inverse_powers(Q):
                ys.append(mask)
        for x in xs:
            for y in ys:
                if x & y or x+y == d or gcd(x, y) < 3:
                    continue
                ta, tb = inverse_powers(P), inverse_powers(Q)
                joint = brute_joint(ta[x % P], len(ta), tb[y % Q], len(tb))
                if joint is None:
                    continue
                for u in powers(gcd(x, y)):
                    if u == 2 or d % u == 0:
                        continue
                    lookup = inverse_powers(u)
                    if -d % u in lookup:
                        orbit = brute_joint(*joint, lookup[-d % u], len(lookup))
                        if orbit:
                            found.add(Shape(P, Q, x//Q, y//P, u, *orbit))
    return found


def layer_ok(s, p, w, value):
    x, y, z = s.slopes
    row = x*value, y*value, z*value+1
    n = sum(row)
    divides_scale = p == s.u and n % w == 0
    return (n//w == sum(t//w for t in row)
            and (not divides_scale or (s.k % w == 0 and s.l % w == 0)))


def check_empty_by_exponents(s, p, levels):
    classes = [s.first]
    stride = s.period
    w = 1
    for _ in range(levels):
        w *= p
        step = pow(2, stride, w)
        multiplier, branches = step, 1
        while multiplier != 1:
            multiplier = multiplier*step % w
            branches += 1
        following = []
        for e in classes:
            value = pow(2, e, w)
            for j in range(branches):
                if layer_ok(s, p, w, value):
                    following.append(e+j*stride)
                value = value*step % w
        classes = following
        stride *= branches
    assert not classes, (s, p, levels)


def no_termination(s):
    p = s.u
    w = 1
    while pow(2, s.period, p*w) == 1:
        w *= p
    v = pow(2, s.first, w)
    x, y, z = s.slopes
    graph = CarryGraph(p, (x,y,z), (x*v//w,y*v//w,(z*v+1)//w))
    assert not graph.accept and len(graph.states) == 229


def direct_finite(s, limit):
    # Direct enumeration of e, with an independent floor test at every layer.
    u = s.u
    w = 1
    while True:
        if pow(2, s.period, w*u) != 1:
            break
        w *= u
    precision = 1
    period = s.period
    while period <= limit-s.first:
        precision *= u
        if precision > w:
            period *= u
    leaves = filtered = 0
    for e in range(s.first or s.period, limit+1, s.period):
        value = pow(2, e, precision)
        mod = u
        while mod <= precision and layer_ok(s, u, mod, value):
            mod *= u
        if mod <= precision:
            continue
        leaves += 1
        if not all(layer_ok(s, p, power*p, pow(2, e, power*p))
                   for p, power in zip(s.bases[1:], (s.P, s.Q))):
            continue
        filtered += 1
        # If a coefficient survives this far, the individual exponent must
        # still fail a true factorial-floor test at one of the forced bases.
        assert any(not layer_ok(s, p, p**h, pow(2, e, p**h))
                   for p in s.bases for h in range(1, 65))
    return leaves, filtered


@lru_cache(None)
def can_finish(row, state):
    graph = CarryGraph(3, row, state)
    return bool(graph.accept)


def direct_companion(s, limit):
    # Enumerate exponents directly and check the normalized companion's
    # digit sums by factorial floors at every retained ternary layer.
    precision, period = 81, s.period
    while period <= limit-s.first:
        precision *= 3
        period *= 3
    leaves = filtered = 0
    row = (8, 15, 24)
    for e in range(s.first or s.period, limit+1, s.period):
        t = pow(2, e, precision)
        v = (s.P*s.Q*t+1)*pow(47, -1, precision) % precision
        parts = tuple(a*v for a in row)
        w = 3
        while w <= precision and sum(parts)//w == sum(a//w for a in parts):
            w *= 3
        if w <= precision or not can_finish(row, tuple(a//precision for a in parts)):
            continue
        leaves += 1
        filtered += all(layer_ok(s, p, power*p, pow(2, e, power*p))
                        for p, power in zip(s.bases[1:], (s.P, s.Q)))
    return leaves, filtered


def check_short_split(s, limit):
    # Determine second-digit tests by factorial floors on literal exponent classes.
    for p, power in zip(s.bases[1:], (s.P, s.Q)):
        w = p*power
        ratio = pow(2, s.period, w)
        if ratio == 1:
            continue
        assert pow(ratio, p, w) == 1
        allowed = [j for j in range(p)
                   if layer_ok(s, p, w, pow(2, s.first+j*s.period, w))]
        if len(allowed) <= 2:
            totals = [0, 0]
            for j in allowed:
                sub = s._replace(first=s.first+j*s.period, period=s.period*p)
                a, b = check_short_split(sub, limit)
                totals[0] += a
                totals[1] += b
            return tuple(totals)
    if (s.P, s.Q) == (379, 797):
        return direct_companion(s, limit)
    return direct_finite(s, limit)


def main():
    for base, z in ((3, Fraction(4, 3)), (5, Fraction(6, 5))):
        assert all(z**s * sum(z**(-j) for j in range(base-s)) <= Fraction(3*base+1, 4)
                   for s in range(base))
    counts, shapes = generate(1000)
    assert counts['pairs'] == 16805 and len(shapes) == len(set(shapes)) == 288
    assert set(shapes) == mask_shapes(1000)
    remaining = []
    local = terminal = 0
    for s in shapes:
        for p in s.bases:
            depth = first_empty_level(s, p, 24, 200000)
            if depth is not None:
                check_empty_by_exponents(s, p, depth)
                local += 1
                break
        else:
            if terminal_possible(s):
                remaining.append(s)
            else:
                no_termination(s)
                terminal += 1
    assert (local, terminal, len(remaining)) == (277, 1, 10)
    assert len({(s.P, s.Q) for s in remaining}) == 9
    print('All 288 classes matched by submasks, 277 local exclusions checked by exponent orbits, '
          'one 229-state terminal exclusion.', flush=True)
    leaves = filtered = 0
    for s in remaining:
        limit = 10**8
        expected = direct_finite(s, limit)
        tables = [predecessor_table(s, w) for w in (s.P, s.Q)]
        assert finite_progression(s, limit, tables) == expected
        # Also exercise the early sparse predecessor-class optimization.
        assert finite_check(s, limit) == check_short_split(s, limit)
        leaves += expected[0]
        filtered += expected[1]
    special, = (s for s in remaining if (s.P, s.Q) == (379, 797))
    assert finite_check(special, 10**12) == check_short_split(special, 10**12)
    ternary = CarryGraph(3, (8, 15, 24), (0, 0, 0))
    assert 2481 % 81 == 51 and ternary.accepts_integer(2481)
    assert leaves > filtered > 0
    print(f'Direct exponent controls through {limit}: {leaves} scale leaves, '
          f'{filtered} passing both predecessor layers.')


if __name__ == '__main__':
    main()
