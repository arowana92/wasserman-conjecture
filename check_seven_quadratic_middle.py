"""Certify the binary orders allowed by the middle digits of the quadratic rows."""

from fractions import Fraction
from math import isqrt, lcm

from check_seven_quadratic_bases import limits
from check_seven_quadratic_suffixes import RESIDUAL, upper_graph
from check_seven_ternary_words import offsets


def shift_prefix(D, f, parity, depth):
    square = 100 * Fraction(8*D, 5) * f * 2**parity
    numerator = square.numerator << (2*depth)
    floor = isqrt(numerator // square.denominator)
    # Outward margin protects the digits when 10 alpha is replaced by Y / 2^u.
    margin = 1 << 400
    assert (floor*margin+1)**2 * square.denominator < numerator*margin**2
    assert numerator*margin**2 < ((floor+1)*margin-1)**2 * square.denominator
    return floor, floor.bit_length()-1-depth


def transitions(D, z, w):
    slopes, states = offsets(D, z, w)
    ids = {s: i for i, s in enumerate(states)}
    upper = [[set() for _ in states] for _ in (0, 1)]
    forward = [[set() for _ in states] for _ in (0, 1)]
    terminal = set()
    for i, (x, y, v) in enumerate(states):
        if v >= 0 and not (x & y or x & v or y & v):
            terminal.add(i)
        for bit in (0, 1):
            a, b, c = (t+bit*s for t, s in zip((x, y, v), slopes))
            for shift in (0, 1):
                if (a & 1) + (b & 1) + ((c-shift) & 1) > 1:
                    continue
                j = ids[a//2, b//2, (c-shift)//2]
                forward[shift][i].add(j)
                if shift == 0:
                    upper[bit][j].add(i)
    return terminal, upper, forward


def order_classes(D, f, parity, terminal, upper, forward):
    shift, t = shift_prefix(D, f, parity, 200)
    viable = set(range(len(forward[0])))
    for k in range(shift.bit_length()):
        transitions = forward[(shift >> k) & 1]
        viable = {j for i in viable for j in transitions[i]}

    prefix_length = 32
    prefix = (f.numerator << prefix_length) // f.denominator
    assert prefix.bit_length() == prefix_length
    fraction = f*2**prefix_length - prefix
    assert min(fraction, 1-fraction) > Fraction(210, 2**(100-prefix_length))
    group = frozenset(terminal)
    for digit in bin(prefix)[2:]:
        group = frozenset(j for i in group for j in upper[int(digit)][i])

    history, allowed = {}, []
    while group not in history:
        history[group] = len(allowed)
        allowed.append(bool(group & viable))
        group = frozenset(j for i in group for digit in (0, 1) for j in upper[digit][i])
    assert history[group] == 0
    period = len(allowed)
    # h = u + parity + 4 - t; j = u + floor(log2 alpha).
    alpha_square = Fraction(8*D, 5) * f * 2**parity
    alpha_order = isqrt(alpha_square.numerator//alpha_square.denominator).bit_length()-1
    assert 4**alpha_order < alpha_square < 4**(alpha_order+1)
    residues = {(prefix_length+k+t-parity-4+alpha_order) % period
                for k, ok in enumerate(allowed) if ok}
    return period, residues


def check():
    results = {}
    union = {}
    for D, _, _, z, w in RESIDUAL:
        terminal, upper, forward = transitions(D, z, w)
        _, graph = upper_graph(D, z, w)
        _, rays = limits(graph)
        cases = []
        for f in rays:
            for parity in (0, 1):
                period, residues = order_classes(D, f, parity, terminal, upper, forward)
                results[D, f, parity] = period, residues
                cases.append((period, residues))
                print(D, f, parity, period, sorted(residues))
        period = lcm(*(p for p, _ in cases))
        union[D] = period, {j for j in range(period) if any(j % p in r for p, r in cases)}
    common = lcm(*(p for p, _ in union.values()))
    eligible = {j for j in range(common) if any(j % p in r for p, r in union.values())}
    print('Eligible orders by cofactor:', [(D,p,len(r)) for D,(p,r) in union.items()])
    print('Union:', len(eligible), 'of', common)
    assert union[4587] == (138, {4, 26, 67, 90})
    assert union[9453] == (68, {15, 23, 40, 64})
    assert union[7383][0] == 1166 and len(union[7383][1]) == 683
    assert (common, len(eligible)) == (2735436, 1698764)
    return results


if __name__ == '__main__':
    check()
