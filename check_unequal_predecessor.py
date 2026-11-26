"""Finite binary gap certificate used in equal_predecessor_gaps.md."""

from argparse import ArgumentParser
from collections import Counter
from math import gcd

from flint import fmpz


def large_pairs():
    """Generate the reduced-denominator possibilities of the two half-integer relations."""
    result = set()
    for d in range(2002, 30001, 4):
        for denominator in range(4, 101, 8):
            divisor, remainder = divmod(2 * d, denominator)
            if remainder:
                continue
            start = d // 4 + 1
            for sign in (-1, 1):
                residue = sign * pow(2, -1, divisor) % divisor
                for e in range(start + (residue - start) % divisor, d, divisor):
                    g = gcd(d, e)
                    if e % 4 != 2 or g > 28 or d >= 3 * e + g:
                        continue
                    if any(2 * e // gcd(2 * e, 2 * d + b) <= 100 for b in (-1, 1)):
                        result.add((d, e))
    return result


def large_pairs_by_slopes():
    """Solve 2 A e = 2 B d + sign*A at both rational slopes."""
    result = set()
    slopes = [(b, a) for a in range(4, 101, 8) for b in range(1, a, 2) if gcd(b, a) == 1]
    inverses = [(b, a) for a in range(4, 101, 8) for b in range(a + 1, 4 * a, 2)
                if gcd(b, a) == 1]
    for b, a in slopes:
        for f, c in inverses:
            det = 2 * (a * c - b * f)
            assert det
            for s in (-1, 1):
                for t in (-1, 1):
                    d, rem = divmod(a * (f * s + c * t), det)
                    if rem or d <= 2000:
                        continue
                    e, rem = divmod(2 * b * d + s * a, 2 * a)
                    g = gcd(d, e)
                    if not rem and d % 4 == e % 4 == 2 and 0 < e < d < 3 * e + g and g <= 28:
                        assert d < 30000
                        result.add((d, e))
    return result


def eligible_pairs(limit):
    for d in range(2, limit + 1):
        for e in range(d // 4 + 1, d):
            g = gcd(d, e)
            if g <= 28 and (d // g) % 2 and (e // g) % 2 and d < 3 * e + g:
                yield d, e


def inspect(pairs):
    counts = Counter()
    survivors = []
    for d, e in pairs:
        g = gcd(d, e)
        common = (1 << g) + 1
        d_power = fmpz(1) << d
        e_power = fmpz(1) << e
        d_mod, e_mod = (d_power + 1) // common, (e_power + 1) // common
        counts['gap pairs'] += 1
        d_cap = (pow((e_power + 1) % d_mod, 2 * d, d_mod) - 1).gcd(d_mod)
        threshold = 27 * d_power**2
        if common**3 * (d_cap * e_mod)**2 < threshold:
            continue
        counts['first order bound'] += 1
        e_cap = (pow((d_power + 1) % e_mod, 2 * e, e_mod) - 1).gcd(e_mod) if e_mod > 1 else 1
        if common**3 * (d_cap * e_cap)**2 >= threshold:
            survivors.append((d, e))
            exact_d = max(gcd(int(d_cap), (1 << j) * (1 + (1 << e)) - 1) for j in range(2 * d))
            exact_e = max(gcd(int(e_cap), (1 << i) * (1 + (1 << d)) - 1) for i in range(2 * e))
            assert common**3 * (exact_d * exact_e)**2 < threshold
            print('Last pair', d, e, 'order bounds', d_cap, e_cap, 'cyclic bounds', exact_d, exact_e)
    return counts, survivors


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--compare-slopes', action='store_true')
    args = parser.parse_args()
    large = large_pairs()
    if args.compare_slopes:
        assert large == large_pairs_by_slopes()
    for name, pairs in [('initial', eligible_pairs(2000)), ('large', sorted(large))]:
        counts, left = inspect(pairs)
        print(name, dict(counts), 'before cyclic bounds', left, flush=True)
