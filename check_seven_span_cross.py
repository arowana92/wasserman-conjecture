"""Direct integer checks on the normalized allocations and the aligned binary carry."""

from itertools import combinations, combinations_with_replacement, product
from math import gcd, lcm

from check_six_bit_rays import normalized_rows
from check_seven_span import binary_scale, pair_at_one, prime_powers


def direct(base, p, positions, powers):
    if min(positions) or max(positions.count(i) for i in positions) >= p:
        return None
    u, v, j, k, l, singleton = positions
    x = powers[u]+powers[v]
    y = powers[j]+powers[k]+powers[l]
    z = powers[singleton]
    m = x+y+z
    g = gcd(x, y)
    for q in (2, p):
        while g % q == 0:
            g //= q
    eligible = gcd(m, lcm(x, y))
    while eligible % 2 == 0:
        eligible //= 2
    if g >= (3 if p == 2 else 7) and 27*m*m <= g**3*eligible**2:
        return positions, m, g, eligible


def bounded(span, base, p):
    powers = [base**i for i in range(span+1)]
    choices = range(span+1)
    triples = combinations if p == 2 else combinations_with_replacement
    seen = set()
    tests = 0
    for pair, trio, k in product(combinations(choices, 2), triples(choices, 3), choices):
        if trio[0] == trio[-1]:
            continue
        result = direct(base, p, pair+trio+(k,), powers)
        if result:
            seen.add(result)
        tests += 1
    assert seen == set(pair_at_one(base, p, span))
    return tests


def selected(span, base, p, gaps):
    d, e, f = gaps
    powers = [base**i for i in range(span+1)]
    expected = set()
    for i in range(span-d+1):
        for j in range(span-f+1):
            for k in (range(span+1) if min(i, j) == 0 else (0,)):
                result = direct(base, p, (i, i+d, j, j+e, j+f, k), powers)
                if result:
                    expected.add(result)
    actual = {v for v in pair_at_one(base, p, span)
              if (v[0][1]-v[0][0], v[0][3]-v[0][2], v[0][4]-v[0][2]) == gaps}
    assert expected == actual and expected
    return len(expected)


def aligned_example():
    positions = (0, 85, 1, 86, 3, 87)
    f = 2**85+1
    m = sum(2**i for i in positions)
    # k+h is the upper singleton position of the shifted triple.
    # Here the smaller position of that triple must be 86 modulo 170.
    shift = 86-3
    assert (1+(2**shift)*(2**3+2**87)) % f == 0
    reduced = f
    for q, w in prime_powers(f).items():
        first = (2**shift)*f//w
        second = 2*first
        layer = q
        while layer <= w:
            if first % layer + second % layer >= layer:
                reduced //= w
                break
            layer *= q
    assert reduced == binary_scale(positions, f) == f//3
    assert not list(normalized_rows(reduced, tuple(prime_powers(m).values())))
    return reduced


if __name__ == '__main__':
    print('Direct allocation controls', sum(bounded(9, b, p) for b, p in ((2, 2), (4, 2), (3, 3), (5, 5))))
    templates = [(2, 2, (8, 7, 23)), (4, 2, (27, 18, 36)), (3, 3, (18, 12, 24))]
    print('Matched positive controls at span 100', sum(selected(100, *args) for args in templates))
    print('Reduced scale for the 85-bit gap', aligned_example())
