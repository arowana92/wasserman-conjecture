"""Count possible orders of magnitude in a window with arbitrary starting index."""
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from math import gcd

from carry_blocks import CarryGraph, PrefixMatcher
from check_carry_intersections import dimension_percent
from check_power_obstruction_cover import DENOMINATORS
from check_prime_power_cofactors import shapes
from check_shifted_directions import possible_scales
from search_endpoint import atoms, smallest_factors


def ceiling(x, y):
    return -(-x // y)


@lru_cache(maxsize=None)
def log_interval(a, b, bits):
    """Outward dyadic enclosure of log(a/b) on [1,2], in units 2^-bits."""
    assert b <= a <= 2*b
    scale = 1 << bits
    num, den = a-b, a+b
    term_lo, term_hi = scale*num//den, ceiling(scale*num,den)
    lower = upper = 0
    for j in range(bits//3+2):
        lower += (2*term_lo)//(2*j+1)
        upper += ceiling(2*term_hi,2*j+1)
        term_lo = (term_lo*num*num)//(den*den)
        term_hi = ceiling(term_hi*num*num,den*den)
    # With t <= 1/3 the uncounted tail is below 9^(-number of terms).
    return lower, upper+1


@lru_cache(maxsize=None)
def log_integer(n,bits):
    e = n.bit_length()-1
    lo,hi = log_interval(n,1<<e,bits)
    two_lo,two_hi = log_interval(2,1,bits)
    return lo+e*two_lo, hi+e*two_hi


def log_ratio_bounds(p,r,bits):
    pl,pu = log_integer(p,bits)
    rl,ru = log_integer(r,bits)
    return Fraction(pl,ru), Fraction(pu,rl)


@lru_cache(maxsize=None)
def rotation_denominator(p,r,H):
    lo,hi = log_ratio_bounds(p,r,2*H.bit_length()+100)
    midpoint = (lo+hi)/2
    n,d = midpoint.numerator,midpoint.denominator
    a0,a1,b0,b1 = 0,1,1,0
    best = None
    while d:
        term,rem = divmod(n,d)
        a0,a1 = a1,term*a1+a0
        b0,b1 = b1,term*b1+b0
        if b1 > H:
            break
        if max(abs(b1*lo-a1),abs(b1*hi-a1))*b1 < 1:
            best = a1,b1
        n,d = d,rem
    assert best is not None and gcd(*best)==1
    return best[1]


def prefix_counts(graph,depth):
    """Numbers of distinct nonzero-leading words with any reachable lower state."""
    read = PrefixMatcher(graph,1)
    counts = {graph.accept:1}
    transitions = {}
    totals = [1]
    for i in range(depth):
        following = Counter()
        for mask,ways in counts.items():
            if mask not in transitions:
                transitions[mask] = [read.previous(mask,d) for d in range(graph.base)]
            for dest in transitions[mask][int(i==0):]:
                if dest:
                    following[dest] += ways
        counts = following
        totals.append(sum(counts.values()))
    return totals


@lru_cache(maxsize=None)
def graph_for(p,slopes):
    return CarryGraph(p,slopes,(0,0,1))


def exponent_floor(n,p):
    e = 0
    while n >= p:
        n //= p
        e += 1
    return e


@lru_cache(maxsize=None)
def forced_configurations():
    spf = smallest_factors(1000)
    found = set()
    for den in DENOMINATORS:
        factors = atoms(den, spf)
        for shape in shapes(den, factors):
            for _, _, bases in possible_scales(den, factors, shape, spf):
                found.add((den, shape[1:4], tuple(sorted(bases))))
    return tuple(sorted(found))


def occupied_bound(p, r, H, left_counts, right_counts):
    """An upper bound on occupied base-p orders in any block of H orders."""
    b = rotation_denominator(p, r, H)
    n0 = exponent_floor(b, p) + 1
    m0 = exponent_floor(b, r) + 1
    ns = range(max(1, n0 - 3), n0 + 4)
    ms = range(max(1, m0 - 3), m0 + 4)
    left = left_counts(max(ns))
    right = right_counts(max(ms))
    best = H
    assert r >= 3 and log_integer(r, 32)[0] >= 1 << 32
    # The integer lower bound on log(r) keeps the visit count outward.
    log_r_lower = max(1, (r.bit_length() - 1) * 69 // 100)
    assert log_integer(r, 32)[0] >= log_r_lower << 32
    for n in ns:
        for m in ms:
            P, R = p**(n - 1), r**(m - 1)
            arcs = left[n] * right[m]
            visits = 3*arcs + (b*arcs*(P + R)) // (P*R*log_r_lower)
            small = exponent_floor(max(P, R), p) + 1
            best = min(best, small + ceiling(H, b)*min(b, visits))
    return best


@lru_cache(maxsize=None)
def dimension(p, slopes):
    return dimension_percent(graph_for(p, slopes))


def forced_pairs():
    pairs = set()
    for den, slopes, bases in forced_configurations():
        dims = sorted((dimension(p, slopes), p) for p in bases)
        assert dims[0][0] + dims[1][0] <= 86
        p, r = sorted((dims[0][1], dims[1][1]))
        assert r <= 19
        pairs.add((den, slopes, p, r))
    return sorted(pairs)


def row_bound(slopes, p, r, H):
    return occupied_bound(p, r, H,
                          lambda n: prefix_counts(graph_for(p, slopes), n),
                          lambda n: prefix_counts(graph_for(r, slopes), n))


def verify():
    from math import comb
    families = forced_pairs()
    assert len(families) == 16
    cubes = ((119, 203, 664), (203, 714, 69))
    for slopes, j, k in ((cubes[0], 287, 439), (cubes[1], 251, 412)):
        for base, s in ((7, j), (17, k)):
            assert graph_for(base, slopes).block_bound(1024)**1000 <= base**(1024*s)

    binary_cut = CarryGraph(2, (6, 1), (0, 1))
    assert binary_cut.block_bound(1024).bit_length() == 569
    for k in (12, 30, 100):
        H = 10**k
        cube = min(H, sum(row_bound(sl, 7, 17, H) for sl in cubes))
        all_shapes = min(H, sum(((p-1).bit_length()+1)*row_bound(sl, p, r, H+2)
                                for _, sl, p, r in families))
        # The nonzero base-seven digits increase weakly in reading order from the top.
        seven = occupied_bound(2, 7, H, lambda n: prefix_counts(binary_cut, n),
                               lambda m: [comb(i+5, 5) for i in range(m+1)])
        print(f"H=10^{k}: cube base-7 orders <= {cube}, all-shape binary orders"
              f" <= {all_shapes}, even factor-seven binary orders <= {seven}", flush=True)
        if k == 100:
            assert cube < 166*10**72 and all_shapes < 308*10**78


if __name__ == '__main__':
    verify()
