"""Check the arithmetic ingredients in the sparse-weight counting bounds."""

from collections import Counter
from itertools import combinations, combinations_with_replacement, product
from math import gcd, lcm, prod

from check_three_factors import valuation
from search_endpoint import atoms, smallest_factors


def local_assignments():
    counts = Counter()
    spf = smallest_factors(210)
    for n in range(3, 211):
        for prime in (2, 3, 5, 7):
            if n % prime > 1:
                continue
            predecessor = n % prime == 1
            base = n-1 if predecessor else n
            w = next(power for q, power in atoms(base, spf) if q == prime)
            m = base//w
            for a in range(1, n//3+1):
                for b in range(a, (n-a)//2+1):
                    row = (a, b, n-a-b)
                    if valuation(n, row, prime):
                        continue
                    counts['rows'] += 1
                    if predecessor:
                        unshifted = [i for i in range(3) if row[i] % w == 0]
                        assert len(unshifted) == 2
                        shifted = next(i for i in range(3) if i not in unshifted)
                        x, y = (row[i]//w for i in unshifted)
                        z = (row[shifted]-1)//w
                        assert m == x+y+z
                    else:
                        assert all(part % w == 0 for part in row)
                    for q, power in atoms(n-1, spf):
                        if q == prime or valuation(n, row, q):
                            continue
                        zeros = [i for i in range(3) if row[i] % power == 0]
                        assert len(zeros) == 2
                        third = next(i for i in range(3) if i not in zeros)
                        assert (row[third]-1) % power == 0
                        if not predecessor:
                            counts['weight gcds'] += 1
                        elif third == shifted:
                            assert gcd(x, y, z) % power == 0
                            counts['unshifted pair'] += 1
                        else:
                            common = x if unshifted[0] in zeros else y
                            assert gcd(m, common, 1+w*z) % power == 0
                            counts['pair including shifted part'] += 1
    assert min(counts.values()) > 0 and len(counts) == 4
    print('Full powers on carry-free rows:', dict(counts))


def period_checks():
    counts = Counter()
    for p in (2, 3, 5):
        pairs = list(combinations(range(6), 2))
        for first, second in product(pairs, repeat=2):
            d = first[1]-first[0]
            x, y = (sum(p**h for h in pair) for pair in (first, second))
            for tail in combinations_with_replacement(range(5), 2):
                z = sum(p**h for h in tail)
                for a in (1, 7, 20, 83):
                    u = gcd(p**a*x, p**a*y, 1+p**a*z)
                    reduced = a % (2*d)
                    reduced_gcd = gcd(p**reduced*x, p**reduced*y, 1+p**reduced*z)
                    while reduced_gcd % p == 0:
                        reduced_gcd //= p
                    assert u == reduced_gcd
                    counts['comparisons'] += 1
                    counts['nontrivial gcd'] += u >= 3
                    counts['reduced to zero'] += (a != reduced and reduced == 0)
    assert min(counts.values()) > 0
    print('Two-pair shift reduction:', dict(counts))


def spacing_checks():
    totals = Counter()
    # c=2, delta=1/2, k=8, eta=1/16 give k*delta-binomial(k,2)*eta > c.
    for p, h, z in product((2, 3), (256, 512, 1024), (1, 7)):
        chosen = (h//2, h//2+1, h//2+3)
        modulus = prod(1+p**a*z for a in chosen)
        assert modulus <= p**(2*h)
        large = {a: value for a in range(1, 2*h+1)
                 if (value := gcd(modulus, 1+p**a*z)) >= p**(h//2)}
        assert set(chosen) <= large.keys()
        for a in large:
            window = [(i, v) for i, v in large.items() if a <= i <= a+h//16]
            assert len(window) < 8
            pair_gcds = []
            for (i, u), (j, v) in combinations(window, 2):
                common = gcd(u, v)
                assert (p**abs(i-j)-1) % common == 0
                pair_gcds.append(common)
            ds = [v for _, v in window]
            assert lcm(*ds)*prod(pair_gcds) >= prod(ds)
            totals['windows with multiple hits'] += len(window) > 1
        totals['shifts'] += 2*h
        totals['large gcds'] += len(large)
    assert min(totals.values()) > 0
    print('Large-divisor spacing:', dict(totals))


def seven_shift_controls():
    shapes = ((2, (0, 1, 2, 6, 7, 8)), (3, (0, 0, 2, 2, 4, 4)),
              (5, (0, 0, 1, 1, 2, 2)), (2, (0, 1, 2, 100, 101, 3000)))
    counts = Counter()
    for p, shape in shapes:
        powers = tuple(p**e for e in shape)
        m, span = sum(powers), max(shape)
        possible = set()
        allocations = 0
        for pair in combinations(range(6), 2):
            remain = tuple(i for i in range(6) if i not in pair)
            x = sum(powers[i] for i in pair)
            for size in (2, 3):
                for chosen in combinations(remain, size):
                    if size == 2 and chosen < pair:
                        continue
                    allocations += 1
                    y = sum(powers[i] for i in chosen)
                    z = m-x-y
                    scale = gcd(x, y)
                    while gcd(scale, p*m) > 1:
                        scale //= gcd(scale, p*m)
                    if scale < 3:
                        continue
                    d = gcd(x, y, z)
                    reduced = m//d
                    assert (2*reduced)**3 >= p**(2*span)
                    k = reduced * gcd(x//d, y//d)
                    shifts = [a for a in range(1, 341)
                              if gcd(k, 1+p**a*z)**3 >= 27*reduced**2]
                    possible.update(shifts)
                    for left, right in zip(shifts, shifts[3:]):
                        assert 9*(right-left) > span
                        assert p**(12*(right-left)) >= 3**12 * reduced**2
                        counts['four-shift windows'] += 1
                    counts['positive shift tests'] += len(shifts)
        assert allocations == 105
        counts['labelled allocations'] += allocations
        ordered = sorted(possible)
        assert all(9*(b-a) > span for a, b in zip(ordered, ordered[315:]))
    assert min(counts.values()) > 0
    print('Explicit seven-unit shift bound:', dict(counts))


def positive_controls():
    count = 0
    for p, k, j in product((2, 3, 5, 7, 11), range(1, 31), range(1, 12)):
        t = p**k
        pair = p**j*(1+t**3)
        triple = 1+t**2+t**4
        assert gcd(pair, triple) == gcd(pair, p**(j+2)*(1+t**3), triple) == t*t-t+1
        count += 1
    print('Exponentially large gcd controls:', count)


if __name__ == '__main__':
    local_assignments()
    period_checks()
    spacing_checks()
    seven_shift_controls()
    positive_controls()
