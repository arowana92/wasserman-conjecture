"""Construct a progression hiding arbitrary depth at any prescribed finite prime list."""

from argparse import ArgumentParser
from math import gcd, isqrt

D = 89 * 293
SLOPES = (9345, 16408, 324)
# At these four scalars, the complete first row avoids the indicated prime.
SCALARS = {2: 2, 7: 108, 89: 9587, 293: 1654}


def prime(p):
    return p >= 2 and all(p % i for i in range(2, isqrt(p) + 1))


def add_congruence(a, period, b, modulus):
    a += period * (((b - a) * pow(period, -1, modulus)) % modulus)
    return a, period * modulus


def first_row(w):
    a, b, c = SLOPES
    return a*w, b*w, c*w + 1


def passes(row, p, depth=None):
    n = sum(row)
    q = p
    j = 0
    while q <= n and (depth is None or j < depth):
        if n // q != sum(x // q for x in row):
            return False
        q *= p
        j += 1
    return True


def progression(primes, extra, witness):
    chosen = sorted(set(primes) | SCALARS.keys())
    assert all(map(prime, chosen)) and extra >= 1
    assert prime(witness) and witness > max(2*D, *chosen)
    w, period = 0, 1
    for p in chosen:
        modulus = p**(extra + 1)
        if p in SCALARS:
            residue = SCALARS[p]
        else:
            v = p * (2 if p == 5 else 1)
            residue = (7*v - 1) * pow(D, -1, modulus) % modulus
        w, period = add_congruence(w, period, residue, modulus)
    modulus = witness**2
    residue = (witness - 1) * pow(D, -1, modulus) % modulus
    return add_congruence(w, period, residue, modulus)


def verify(primes, extra, witness):
    w, period = progression(primes, extra, witness)
    assert gcd(w//2, period//2) == 1 and w % 4 == 2
    for scalar in (w, w + period):
        n = D*scalar + 1
        a = first_row(scalar)
        v, rem = divmod(n, 7)
        assert rem == 0
        c = (v, 3*v, 3*v)
        assert sum(a) == sum(c) == n
        assert gcd(scalar, D) == 1 and gcd(v, 7) == 1
        for p in set(primes) | SCALARS.keys():
            assert n*(n-1) % p == 0 and n*(n-1) % (p*p) != 0
            assert passes(a if p in SCALARS else c, p, extra+1)
        assert n % (witness**2) == witness
        assert not passes(a, witness, 1) and passes(c, witness, 1)
        assert not passes(c, witness, 2)
    return w, period


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--prime-limit', type=int, default=101)
    parser.add_argument('--extra', type=int, default=16)
    parser.add_argument('--print-progression', action='store_true')
    args = parser.parse_args()
    chosen = [p for p in range(2, args.prime_limit+1) if prime(p)]
    witness = max(2*D, args.prime_limit) + 1
    while not prime(witness):
        witness += 1
    for p, w in SCALARS.items():
        assert passes(first_row(w), p)
    w, modulus = verify(chosen, args.extra, witness)
    print('chosen primes including template:', len(set(chosen) | SCALARS.keys()))
    print('additional digits:', args.extra, 'common witness:', witness)
    print('progression modulus bit length:', modulus.bit_length())
    if args.print_progression:
        print(w, modulus)
