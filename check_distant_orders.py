"""Certify intervals at huge binary exponents using only the leading digits."""
from argparse import ArgumentParser
from functools import lru_cache

from carry_blocks import CarryGraph, PrefixMatcher
from check_logarithmic_intervals import (
    ceiling, forced_configurations, log_ratio_bounds,
)

SCALE_BITS = 512
SCALE = 1 << SCALE_BITS
START = 10**100
LENGTH = 16384
MARGIN_BITS = 440


def product(x, y, bits):
    """Multiply intervals [lo,hi] times a power of two, rounding outward."""
    lo, hi, e = x
    left, right, f = y
    lower, upper = lo*left, hi*right
    chop = max(0, upper.bit_length() - bits)
    return lower >> chop, ceiling(upper, 1 << chop), e+f+chop


def power_interval(p, n, bits):
    result = (1, 1, 0)
    factor = (p, p, 0)
    while n:
        if n & 1:
            result = product(result, factor, bits)
        n //= 2
        if n:
            factor = product(factor, factor, bits)
    return result


def root_bound(p, j):
    """Enclose SCALE*p**e/2**j with p**e >= 2**(j+1)."""
    if p == 2:
        return 2*SCALE, 2*SCALE, j+1
    l, h = log_ratio_bounds(2, p, j.bit_length()+150)
    target = (j+1)*l
    e = target.numerator // target.denominator + 1
    assert e-1 < target <= (j+1)*h < e
    a, b, f = power_interval(p, e, j.bit_length()+SCALE_BITS+64)
    shift = f-j+SCALE_BITS
    if shift >= 0:
        lo, hi = a << shift, b << shift
    else:
        lo, hi = a >> -shift, ceiling(b, 1 << -shift)
    assert SCALE <= lo <= hi <= 3*p*SCALE and hi-lo <= 2
    return lo, hi, e


def advance(p, current):
    lo, hi, e = current
    if p == 2:
        return 2*SCALE, 2*SCALE, e+1
    # At the new index, the old normalized power has halved. Multiply
    # it by p when needed so every possible input still fits in e digits.
    if lo < 4*SCALE:
        lo, hi, e = lo*p, hi*p, e+1
    lo, hi = lo//2, ceiling(hi, 2)
    assert 2*SCALE <= lo <= hi <= 3*p*SCALE
    return lo, hi, e


def widened_root(p, bounds):
    """Room for every relative phase perturbation of size 2**(-MARGIN_BITS)."""
    lo, hi, e = bounds
    if p == 2:
        return lo, hi, e
    # One extra leading digit guarantees that every perturbed root still
    # covers the entire binary bin, including a near equality of powers.
    lo, hi = p*lo, p*hi
    radius = ceiling(hi, 1 << MARGIN_BITS)
    assert lo-radius >= 2*SCALE
    return lo-radius, hi+radius, e+1


def cover(bounds, matchers):
    """True only if no simultaneous carry-free input can belong to this bin."""
    cells = tuple((0, 1, m.graph.accept, 0) for m in matchers)
    stack = [(SCALE, 2*SCALE, cells)]
    nodes = 0
    while stack:
        lo, hi, cells = stack.pop()
        nodes += 1
        i = 0
        for k in range(1, len(cells)):
            if bounds[k][1]*cells[i][1] > bounds[i][1]*cells[k][1]:
                i = k
        prefix, denom, mask, depth = cells[i]
        lower, upper, e = bounds[i]
        matcher = matchers[i]
        p = matcher.graph.base
        denom *= p
        if depth >= e or upper < denom:
            return False, nodes
        first = max(0, lo*denom//upper - prefix*p - 1)
        last = min(p-1, hi*denom//lower - prefix*p)
        for digit in range(first, last+1):
            prev = matcher.previous(mask, digit)
            if not prev:
                continue
            value = prefix*p+digit
            left = max(lo, lower*value//denom)
            right = min(hi, ceiling(upper*(value+1), denom))
            if left <= right:
                child = (value, denom, prev, depth+1)
                stack.append((left, right, cells[:i]+(child,)+cells[i+1:]))
    return True, nodes


@lru_cache(maxsize=None)
def matcher(p, slopes, offsets):
    return PrefixMatcher(CarryGraph(p, slopes, offsets), 1)


def small_controls():
    from math import comb
    for p in (2, 3, 7, 17, 29):
        root = root_bound(p, 0)
        for j in range(1025):
            lo, hi, e = root
            assert (lo << j) <= SCALE*p**e <= (hi << j)
            assert p**e >= 2**(j+1)
            root = advance(p, root)

    passed = excluded = 0
    rows = ((2, (1, 2, 4), (0, 0, 0)), (3, (1, 1), (0, 1)))
    for j in range(8):
        free = False
        for t in range(1 << j, 1 << (j+1)):
            for p, sl, off in rows:
                entries = [a*t+b for a,b in zip(sl,off)]
                total = sum(entries)
                coeff = 1
                for x in entries:
                    coeff *= comb(total, x)
                    total -= x
                if coeff % p == 0:
                    break
            else:
                free = True
                passed += 1
        empty, _ = cover([root_bound(p, j) for p, _, _ in rows],
                         [matcher(*row) for row in rows])
        assert not (empty and free)
        excluded += empty
    assert passed > 0 and excluded > 0
    print(f"Small controls: {passed} direct positive inputs and {excluded} certified empty bins.")


def verify(start, count, include_cuts):
    families = {(sl, (0, 0, 1), ps) for _, sl, ps in forced_configurations()}
    if include_cuts:
        # Both parities in the earlier five/seven predecessor reductions.
        for den, k, u in ((5,3,3), (7,3,3), (7,5,5), (7,6,3)):
            families.add(((k, den-k), (0, 1), (2, u, den)))
    groups = [(ps, [matcher(p, sl, off) for p in ps])
              for sl, off, ps in sorted(families)]
    primes = sorted({p for ps, _ in groups for p in ps})
    print("Bases:", primes, flush=True)
    total = worst = 0
    for step in range(count):
        j = start+step
        if step % 4096 == 0:
            bounds = {p: root_bound(p, j) for p in primes}
        else:
            bounds = {p: advance(p, bounds[p]) for p in primes}
        padded = {p: widened_root(p, bounds[p]) for p in primes}
        for ps, matchers in groups:
            empty, nodes = cover([padded[p] for p in ps], matchers)
            assert empty, (j, ps, matchers[0].graph.slopes)
            total += nodes
            worst = max(worst, nodes)
    print(f"Certified {count} binary bins starting at {start}, {len(groups)} cases each,"
          f" with {total} prefix nodes, maximum {worst}.", flush=True)


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--start', type=int, default=START)
    parser.add_argument('--count', type=int, default=LENGTH)
    parser.add_argument('--include-cuts', action='store_true')
    args = parser.parse_args()
    small_controls()
    verify(args.start, args.count, args.include_cuts)
