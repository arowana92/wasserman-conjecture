"""Audit the two-level families using prime-power partitions at the odd part of X+1."""

from argparse import ArgumentParser
from array import array
from collections import Counter
from itertools import combinations
from math import gcd, isqrt

from check_six_bit_rays import normalized_rows


def factor_sieve(limit):
    # At an odd composite store a prime divisor <= sqrt(limit).
    sieve = array('H', [0]) * (limit//2 + 1)
    for p in range(3, isqrt(limit) + 1, 2):
        if not sieve[p//2]:
            start = p*p//2
            count = (len(sieve)-1-start)//p + 1
            sieve[start::p] = array('H', [p]) * count
    return sieve


def odd_atoms(n, sieve):
    assert n % 2
    result = {}
    while n > 1:
        p = sieve[n//2] or n
        w = 1
        while n % p == 0:
            n //= p
            w *= p
        result[p] = w
    return result


def floor_carry(pair, q, largest):
    a, b = pair
    w = q
    while w <= largest:
        if (a+b)//w > a//w + b//w:
            return True
        w *= q
    return False


def audit(x, delta, sieve, counts):
    counts['profiles'] += 1
    f = x+1
    j = f//(f & -f)
    m = 2*f + x**delta
    if 27*m*m > j**3:
        return
    shifted = odd_atoms(m, sieve)
    if j <= 2*max(shifted.values()):
        return
    # Choose the least positive exponent of the required parity.
    power = x**(1+delta)
    n = m*power+1
    assert n % f == 0
    a = power*f
    first = odd_atoms(j, sieve)
    failed = [w for q, w in first.items() if floor_carry((a, a), q, w*w)]
    for row in normalized_rows(j, tuple(shifted.values())):
        if any(any(r % w for r in row) for w in failed):
            continue
        counts['rows'] += 1
        if any(gcd((r+s)//gcd(r,s), 2*m) > 1 for r,s in combinations(row,2)):
            counts['pair divisor'] += 1
            continue
        actual = [r*(n//j) for r in row]
        if any(floor_carry([z for z in actual if z%w == 0], q, w*w)
               for q, w in shifted.items()):
            counts['predecessor carry'] += 1
            continue
        # All odd x have x^64 = 1 modulo 64. Test every exponent of this parity.
        assert pow(x, 64, 64) == 1
        for e in range(1+delta, 65+delta, 2):
            scale = ((1+m*pow(x, e, 64))*pow(j, -1, 64)) % 64
            parts = [scale*r for r in row]
            if not any(sum(parts)//w > sum(z//w for z in parts)
                       for w in (2, 4, 8, 16, 32, 64)):
                raise ArithmeticError(('unresolved binary class', x, delta, row, e))
        counts['binary orbit'] += 1


def check(limit):
    sieve = factor_sieve(3*limit+3)
    counts = Counter()
    for x in range(5, limit+1, 2):
        if x % 3 == 0:
            continue
        if sieve[x//2] and len(odd_atoms(x, sieve)) != 1:
            continue
        for delta in (0, 1):
            audit(x, delta, sieve, counts)
    print(dict(sorted(counts.items())))


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('limit', type=int, nargs='?', default=5_000_000)
    check(parser.parse_args().limit)
