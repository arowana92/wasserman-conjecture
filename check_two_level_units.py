"""Certify all outer exponents at prime-power gaps through five million."""

from argparse import ArgumentParser
from collections import Counter
from itertools import combinations
from math import gcd, prod

from check_close_power_gaps import pair_carries
from check_six_digit_gaps import complement_by_parts, factors
from search_even_endpoint import prime_sieve


def first_divisor(x, delta):
    factor = x + 1
    sign = 2*delta - 1
    allowed = []
    for q, w in factors(factor).items():
        digit = (sign * (factor//w)) % w
        if q != 2 and not pair_carries(digit, digit, q):
            allowed.append(w)
    return prod(allowed)


def binary_orbit(x, delta, h, m, row):
    """Does every scale in the indicated exponent parity force a low carry?"""
    for bits in range(1, 13):
        modulus = 1 << bits
        start = pow(x, 1-delta, modulus)
        step = (x*x) % modulus
        inverse = pow(h, -1, modulus)
        residue = start
        while True:
            scale = ((m*residue+1)*inverse) % modulus
            a, b, c = ((scale*r) % modulus for r in row)
            if (a | b | c) == a+b+c:
                break
            residue = (residue*step) % modulus
            if residue == start:
                return bits
    return 0


def certify(x, delta, counts):
    counts['profiles'] += 1
    m = (2+delta)*(x+1) + 1-2*delta
    h = first_divisor(x, delta)
    if h**3 < 27*m*m:
        return
    atoms = factors(m)
    if 2*max(atoms.values()) + 1 > h:
        return
    counts['CRT searches'] += 1
    for row in complement_by_parts(h, tuple(atoms.values())):
        counts['rows'] += 1
        if any((a+b)//gcd(a, b) % q == 0 for a, b in combinations(row, 2)
               for q in (2, *atoms)):
            counts['pair divisor'] += 1
            continue
        for q, w in atoms.items():
            tails = [(r//w)*pow(h, -1, w) % w for r in row if r % w == 0]
            if pair_carries(*tails, q):
                counts['predecessor carry'] += 1
                break
        else:
            bits = binary_orbit(x, delta, h, m, row)
            if not bits:
                raise ArithmeticError(('unresolved exponent family', x, delta, h, row))
            counts[f'binary orbit through {bits} bits'] += 1


def check(limit):
    counts = Counter()
    for p in prime_sieve(limit):
        if p == 3:
            continue
        x = p
        while x <= limit:
            for delta in (0, 1):
                certify(x, delta, counts)
            x *= p
    print(dict(sorted(counts.items())))


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('limit', type=int, nargs='?', default=5_000_000)
    check(parser.parse_args().limit)
