"""Bound separations of two exceptional exponents, with no limit on their size."""
from argparse import ArgumentParser
from bisect import bisect_left, bisect_right
from itertools import combinations_with_replacement
from math import gcd
from typing import NamedTuple

from carry_blocks import CarryGraph
from check_binary_exponent_patterns import SMALL
from check_binary_cofactor_rays import live_edges
from check_binary_exponent_progressions import lift_progression, published_rays


class Template(NamedTuple):
    d: int
    p: int
    residue: int
    period: int
    edges: dict


def initial_template(ray):
    row, p, _ = ray
    edges = live_edges(CarryGraph(p, row.slopes, (0,0,1)))
    tight = lift_progression(row, p, edges)
    return Template(row.d, p, tight.a, tight.h, edges)


def carry_state(template, scalar, modulus):
    state = 0
    while modulus > 1:
        scalar, digit = divmod(scalar, template.p)
        state = template.edges[state][digit]
        modulus //= template.p
    return state


def lifted(rows, template, width, period, precision):
    """Read one more scalar digit and maintain exponent classes in sorted order."""
    p = template.p
    multiplier = pow(2, period, precision)
    slope = (multiplier-1)//width % p
    assert slope and multiplier % width == 1
    units = [pow(multiplier, i, precision) for i in range(p)]
    new_gcd = gcd(p*period, template.period)
    first = template.residue % new_gcd
    bins = [[] for _ in range(p)]
    for exponent, state, scalar in rows:
        # (scalar * multiplier^i) gains i * scalar * slope at this digit.
        digit = scalar//width % p
        inverse = pow(scalar*slope % p, -1, p)
        for new_digit, child in template.edges[state].items():
            i = (new_digit-digit)*inverse % p
            e = exponent+i*period
            if e % new_gcd == first:
                q = scalar*units[i] % precision if i else scalar
                bins[i].append((e,child,q))
    return [node for group in bins for node in group]


def near_self(rows, period, gap):
    if len(rows) < 2:
        return []
    n = len(rows)
    kept = []
    before = rows[0][0]+period-rows[-1][0]
    for i, node in enumerate(rows):
        after = rows[i+1][0]-node[0] if i+1 < n else rows[0][0]+period-node[0]
        if min(before,after) <= gap:
            kept.append(node)
        before = after
    return kept


def near_other(rows, others, period, gap):
    if not others:
        return []
    values = [node[0] for node in others]
    kept = []
    for node in rows:
        e = node[0]
        left = bisect_left(values,e)-1
        right = bisect_right(values,e)
        before = e-values[left] if left >= 0 else e+period-values[-1]
        after = values[right]-e if right < len(values) else values[0]+period-e
        if min(before,after) <= gap:
            kept.append(node)
    return kept


def certify(left, right, gap, max_layers=160):
    """Stop once no possible two distinct exponents at this gap remain."""
    p = left.p
    assert p == right.p
    period = gcd(left.period, right.period)
    width = p
    assert pow(2,period,p) == 1
    while pow(2,period,width*p) == 1:
        width *= p
    precision = width*p**max_layers

    def root(template):
        e = template.residue % period
        q = pow(2,e,precision)
        return [(e, carry_state(template, q % width, width), q)]

    same = left is right
    a = root(left)
    b = [] if same else root(right)
    peak = 0
    for layers in range(max_layers):
        if period > gap:
            if same:
                a = near_self(a, period, gap)
            else:
                a,b = near_other(a,b,period,gap),near_other(b,a,period,gap)
        if not a or (not same and not b):
            return layers, peak
        peak = max(peak,len(a)+len(b))
        a = lifted(a,left,width,period,precision)
        if not same:
            b = lifted(b,right,width,period,precision)
        width *= p
        period *= p
    raise ArithmeticError(("gap proof did not terminate",left.d,right.d,gap))


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=int, choices=(3,5,7,11))
    args = parser.parse_args()
    active = {row.d for row in SMALL}
    templates = [initial_template(ray) for ray in published_rays()
                 if ray.row.d in active and (args.base is None or ray.prime == args.base)]
    expected = {(r.d,r.a,r.h) for r in SMALL}
    actual = {(r.d,r.residue,r.period) for r in templates}
    assert actual <= expected and (args.base is not None or actual == expected)
    for left, right in combinations_with_replacement(templates,2):
        if left.p != right.p:
            continue
        gap = 10**48 if left is right and left.d == 12031 else 10**15
        layers, peak = certify(left, right, gap)
        print(f"{left.p:2} {left.d} {right.d} {left.residue} {right.residue} "
              f"10^{len(str(gap))-1} {layers} {peak}", flush=True)


if __name__ == "__main__":
    main()
