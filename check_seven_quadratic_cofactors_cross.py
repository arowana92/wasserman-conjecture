"""Independently enumerate candidate parts and replay the cofactor certificates."""

from argparse import ArgumentParser
from collections import defaultdict
from fractions import Fraction
from itertools import combinations, product
from math import gcd, isqrt


def candidates(D):
    by_group = defaultdict(list)
    for z in range(2, D//2 + 1, 2):
        g = gcd(z, D)
        a = D//g
        if 1 < a < D and gcd(a, g) == 1:
            by_group[a].append(z)
    for a, b in combinations(sorted(by_group), 2):
        if gcd(a, b) != 1:
            continue
        for z, w in product(by_group[a], by_group[b]):
            if 2*(z+w) >= D or ((z//gcd(z,w))%2 and (w//gcd(z,w))%2):
                continue
            k = next(i for i in range(1,6) if (i - D*z)%5 == 0)
            l = next(i for i in range(1,6) if (i - D*w)%5 == 0)
            if k + l <= 4:
                yield (D, a, b, z, w)


def no_low_input(key, r, depth):
    D, a, b, z, w = key
    assert r == 2 or (D%r == 0 and all(r%d for d in range(2, isqrt(r)+1)))
    prefixes = [0]
    modulus = 1
    for _ in range(depth):
        last = modulus
        modulus *= r
        inv = pow(D if r == 2 else 10, -1, modulus)
        fixed_a, fixed_b = gcd(a, modulus), gcd(b, modulus)
        passed = []
        for digit in range(r):
            for prefix in prefixes:
                v = digit*last + prefix
                if r == 2:
                    P, R = v, (10*v+7)*inv % modulus
                    if not P%2:
                        continue
                else:
                    R, P = v, (D*v-7)*inv % modulus
                if (2*R*z+5)%fixed_a or (2*R*w+5)%fixed_b:
                    continue
                T = 2*P + 1
                N = T*(5*P+1)
                x, y = T*R*z, T*R*w
                # Legendre's formula at this layer; any integer representatives
                # of the three part residues give the same floor difference.
                if N//modulus == x//modulus + y//modulus + (N-x-y)//modulus:
                    passed.append(v)
        if not passed:
            return True
        prefixes = passed
    return False


def add_without_carry(values, r):
    while any(values):
        if sum(v%r for v in values) >= r:
            return False
        values = [v//r for v in values]
    return True


def no_high_input(key, r, depth):
    D, _, _, z, w = key
    assert r == 2 or (r == 3 and D%3 == 0)
    H = 10**9 + 7
    slopes = (2*z, 2*w, D-2*z-2*w)
    cells = [0]
    for pos in range(depth):
        passed = []
        for previous in cells:
            for digit in range(int(pos == 0), r):
                j = r*previous + digit
                # j <= normalized scale < j+1. Only the last part has error.
                ranges = []
                for i, a in enumerate(slopes):
                    low = a*H - (5*D if i == 2 else 0)
                    start = low*j//(D*H)
                    stop = (a*(j+1)-1)//D
                    ranges.append(range(start, stop+1))
                for row in product(*ranges):
                    if not add_without_carry(row, r):
                        continue
                    lower = [Fraction(j)]
                    upper = [Fraction(j+1)]
                    for i, (v, a) in enumerate(zip(row, slopes)):
                        lower.append(Fraction(D*v, a))
                        upper.append(Fraction(D*H*(v+1), a*H-(5*D if i == 2 else 0)))
                    if max(lower) < min(upper):
                        passed.append(j)
                        break
        if not passed:
            return True
        cells = passed
    return False


def verify(bound, path):
    supplied = {}
    with open(path, encoding='ascii') as data:
        for line in data:
            D, a, b, z, w, reason, prime, depth = line.split()
            key = tuple(map(int, (D,a,b,z,w)))
            assert key not in supplied
            supplied[key] = reason, int(prime), int(depth)
    cases = {key for D in range(1, bound+1) if D%10 in (3,7) for key in candidates(D)}
    assert cases == supplied.keys(), (len(cases), len(supplied))
    remaining = []
    for key in sorted(cases):
        reason, r, depth = supplied[key]
        if reason == 'low':
            assert no_low_input(key, r, depth), (key, r, depth)
        elif reason == 'high':
            assert no_high_input(key, r, depth), (key, r, depth)
        else:
            assert (reason,r,depth) == ('open',0,0)
            remaining.append(key)
    print(f'Independently reconstructed {len(cases)} shapes and validated {len(cases)-len(remaining)} exclusions.')
    print('Remaining:', remaining)
    expected = [(4587,11,417,1668,154), (7383,23,107,642,276), (9453,23,411,822,92)]
    assert remaining == [row for row in expected if row[0] <= bound]


if __name__ == '__main__':
    cli = ArgumentParser(description=__doc__)
    cli.add_argument('bound', type=int)
    cli.add_argument('certificate')
    args = cli.parse_args()
    verify(args.bound, args.certificate)
