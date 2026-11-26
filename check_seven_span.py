"""Certify all seven-unit profiles through a given normalized span, at every shift."""

from argparse import ArgumentParser
from collections import Counter
from functools import lru_cache
from itertools import combinations
from math import gcd, prod

from flint import fmpz

from check_close_power_gaps import pair_carries
from check_six_digit_gaps import complement_by_parts


def odd_part(n):
    return n//(n & -n)


def coprime_part(n, excluded):
    shared = gcd(n, excluded)
    while shared > 1:
        n //= shared
        shared = gcd(n, shared)
    return n


@lru_cache(maxsize=4096)
def prime_powers(n):
    factors = {int(q): int(q**e) for q, e in fmpz(n).factor()}
    assert all(fmpz(q).is_prime() for q in factors) and prod(factors.values()) == n
    return factors


def companion(m, g, eligible, counts):
    assert gcd(m, g) == 1
    forced = {q: w for q, w in prime_powers(m).items() if eligible % w}
    t = prod(forced.values())
    if 27*t*t > g**3 or any(2*w+1 > g for w in forced.values()):
        return
    counts['CRT searches'] += 1
    for row in complement_by_parts(g, tuple(forced.values())):
        counts['CRT triples'] += 1
        if any((i+j)//gcd(i, j) % q == 0 for i, j in combinations(row, 2) for q in forced):
            counts['excluded by a pair divisor'] += 1
            continue
        for q, w in forced.items():
            digits = [((i//w)*pow(g, -1, w)) % w for i in row if i % w == 0]
            if pair_carries(*digits, q):
                counts['excluded by a carry'] += 1
                break
        else:
            raise ArithmeticError(('surviving necessary conditions', m, g, eligible, row))


def caps(modulus, value, powers, span):
    # Bound gcd(modulus, 1+B^r*value) after clearing a unit for negative r.
    return [gcd(modulus, powers[r]*value+1) if r >= 0 else gcd(modulus, value+powers[-r])
            for r in range(-span, span+1)]


def pair_at_one(base, p, span):
    """Yield normalized allocations and exact size-test inputs for the shifted pair."""
    powers = [base**i for i in range(2*span+1)]
    min_scale = 3 if p == 2 else 7
    for d in range(1, span+1):
        f_pair = 1+powers[d]
        for f in range(1, span+1):
            for e in range(1 if p == 2 else 0, f if p == 2 else f+1):
                f_triple = 1+powers[e]+powers[f]
                g = odd_part(gcd(f_pair, f_triple))
                if g < min_scale:
                    continue
                first = coprime_part(odd_part(f_pair), g)
                second = coprime_part(odd_part(f_triple), p*g)
                lower = 27*powers[max(d, f)]**2
                cube = g**3
                if cube*(first*second)**2 < lower:
                    continue
                a_caps = caps(first, f_triple, powers, span)
                a_max = max(a_caps)
                if cube*(a_max*second)**2 < lower:
                    continue
                c_caps = caps(second, f_pair, powers, span)
                if cube*(a_max*max(c_caps))**2 < lower:
                    continue
                for r, c in enumerate(c_caps, -span):
                    if cube*(a_max*c)**2 < lower:
                        continue
                    for t, a in enumerate(a_caps, -span):
                        if cube*(a*c)**2 < lower:
                            continue
                        k = -min(0, r, t)
                        i, j = k+r, k+t
                        positions = (i, i+d, j, j+e, j+f, k)
                        if max(positions) > span or any(positions.count(n) >= p for n in positions):
                            continue
                        m = powers[i]*f_pair+powers[j]*f_triple+powers[k]
                        if cube*(a*c)**2 >= 27*m*m:
                            yield positions, m, g, a*c


def binary_triples(span):
    """All normalized inner shapes after the two aligned gap alternatives."""
    for n in range(1, span+1):
        f = (1 << n)+1
        for h in (n-1, n+1):
            if not 1 <= h <= span:
                continue
            for i, j in combinations(range(span-n+1), 2):
                for k in range(span-h+1):
                    positions = (i, i+n, j, j+n, k, k+h)
                    if min(positions) != 0 or len(set(positions)) != 6:
                        continue
                    m = sum(1 << h for h in positions)
                    if 27*m*m <= f**3:
                        yield positions, m, f


def binary_scale(positions, f):
    i, _, j, _, k, last = positions
    z = (1 << k)+(1 << last)
    scale = 1
    for q, w in prime_powers(f).items():
        multiplier = (-(f//w)*pow(z, -1, w)) % w
        if not pair_carries(((1 << i)*multiplier) % w, ((1 << j)*multiplier) % w, q):
            scale *= w
    return scale


def check(span):
    counts = Counter()
    for p in range(2, 64):
        if not fmpz(p).is_prime():
            continue
        base = p
        while base < 64:
            for positions, m, g, eligible in pair_at_one(base, p, span):
                counts['shifted pair size survivors'] += 1
                companion(m, g, eligible, counts)
            counts['effective bases'] += 1
            base *= p
    print('Shifted pair', dict(counts), flush=True)
    counts = Counter()
    for positions, m, f in binary_triples(span):
        counts['aligned size survivors'] += 1
        g = binary_scale(positions, f)
        if 27*m*m <= g**3:
            counts['after prime-power carry test'] += 1
            companion(m, g, 1, counts)
    print('Shifted triple', dict(counts), flush=True)


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('span', type=int, nargs='?', default=100)
    args = parser.parse_args()
    if args.span < 1:
        parser.error('use a positive span')
    check(args.span)
