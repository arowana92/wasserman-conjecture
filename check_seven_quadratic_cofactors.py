"""Certify all prime-power inputs in a factored quadratic at fixed cofactors."""

from argparse import ArgumentParser
from collections import Counter
from fractions import Fraction
from itertools import combinations, product
from math import gcd

from check_endpoint_pair_bound import prime_powers


def shapes(D):
    bins = [1]
    for _, power in prime_powers(D):
        bins += [v * power for v in bins]
    for a, b in combinations(sorted(bins), 2):
        if a <= 1 or gcd(a, b) != 1:
            continue
        for m in range(2, (a + 1)//2, 2):
            if gcd(m, a) > 1:
                continue
            z = D//a * m
            for n in range(2, (b + 1)//2, 2):
                w = D//b * n
                if 2*(z + w) >= D:
                    break
                if gcd(n, b) > 1 or (z & -z) == (w & -w):
                    continue
                if (D*z - 1)%5 + (D*w - 1)%5 <= 2:
                    yield a, b, z, w


def low_digits(D, shape, r, depth, cap):
    a, b, z, w = shape
    # R = seed + 10v and P = start + Dv parametrize the integral line.
    seed = 7 * pow(D, -1, 10) % 10
    start = (D*seed - 7)//10
    active, place = [0], 1
    for length in range(1, depth + 1):
        layer = r * place
        common_a, common_b = gcd(layer, a), gcd(layer, b)
        kept = []
        for index in active:
            for v in range(index, layer, place):
                R, P = seed + 10*v, start + D*v
                if r == 2 and P%2 == 0:
                    continue
                T = 2*P + 1
                total = T*(5*P + 1) % layer
                first, second = T*R*z % layer, T*R*w % layer
                if (first - 1)%common_a or (second - 1)%common_b:
                    continue
                if first + second <= total:
                    kept.append(v)
                    if len(kept) > cap:
                        return 0
        if not kept:
            return length
        active, place = kept, layer
    return 0


def high_digits(D, shape, r, depth=40, cap=128):
    _, _, z, w = shape
    slopes = (Fraction(2*z, D), Fraction(2*w, D), 1 - Fraction(2*(z+w), D))
    lower = (*slopes[:2], slopes[2] - Fraction(5, 10**9+7))
    assert min(lower) > 0
    digits = [x for x in product(range(r), repeat=3) if sum(x) < r]
    active = [(0, 0, 0)]
    place = 1
    for length in range(1, depth + 1):
        kept = []
        for previous in active:
            for column in digits:
                row = tuple(r*x + y for x, y in zip(previous, column))
                # The normalized scale times 'place' lies in [place,r*place).
                left = max(Fraction(place), *(v/h for v, h in zip(row, slopes)))
                right = min(Fraction(r*place), *((v+1)/h for v, h in zip(row, lower)))
                if left < right:
                    kept.append(row)
                    if len(kept) > cap:
                        return 0
        if not kept:
            return length
        active, place = kept, r*place
    return 0


def exclude(D, shape):
    for r in [2] + [q for q, _ in prime_powers(D)]:
        depth = low_digits(D, shape, r, 96, 256)
        if depth:
            return 'low', r, depth
    for r in (2, 3):
        if r == 3 and D%3:
            continue
        depth = high_digits(D, shape, r)
        if depth:
            return 'high', r, depth
    depth = low_digits(D, shape, 2, 800, 5000)
    if depth:
        return 'low', 2, depth
    return None


def survey(limit, certificate=None):
    used = Counter()
    unresolved = []
    maximum = {}
    for D in range(3, limit+1, 2):
        if D%10 not in (3, 7):
            continue
        for shape in shapes(D):
            result = exclude(D, shape)
            if result:
                end, r, depth = result
                used[end, r] += 1
                maximum[end] = max(maximum.get(end, 0), depth)
            else:
                unresolved.append((D, *shape))
            if certificate is not None:
                print(D, *shape, *(result or ('open', 0, 0)), file=certificate)
    print(f'D through {limit}: {sum(used.values())} exclusions; {len(unresolved)} residual shapes.')
    print('Exclusions by digit end:', {end: sum(v for (side,r),v in used.items() if side == end)
                                       for end in maximum})
    print('Maximum depths:', maximum)
    print('Residual (D, a, b, z, w):', unresolved)
    return unresolved


if __name__ == '__main__':
    cli = ArgumentParser(description=__doc__)
    cli.add_argument('bound', type=int, nargs='?', default=10000, metavar='D')
    cli.add_argument('--certificate', help='write one replayable reason per candidate shape')
    args = cli.parse_args()
    if not 1 <= args.bound <= 10000:
        cli.error('D must be in [1,10000]')
    if args.certificate:
        with open(args.certificate, 'w', encoding='ascii') as output:
            unresolved = survey(args.bound, output)
    else:
        unresolved = survey(args.bound)
    expected = [(4587, 11, 417, 1668, 154), (7383, 23, 107, 642, 276),
                (9453, 23, 411, 822, 92)]
    assert unresolved == [row for row in expected if row[0] <= args.bound]
