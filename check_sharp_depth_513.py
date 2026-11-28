"""Verify the 490-digit weight and its sharp four-prime certificate."""

import gzip
import json
from pathlib import Path

from check_all_prime_depth_witness import profile, prove_chain


def verify():
    h = 512
    u, d = 2**(h+2), 3**(h+1)
    q = d + (-pow(d, -1, u) - d) % u + 233992*u
    n = 1 + d*q
    v, remainder = divmod(n, u)
    r = n - 273
    assert remainder == 0 and q > d and v > u and len(str(n)) == 490
    assert len({2, 3, q, v, r}) == 5
    path = Path(__file__).with_name('depth_513_certificates.json.gz')
    with gzip.open(path, 'rt') as stream:
        certificates = json.load(stream)
    for name, prime in (('q', q), ('v', v), ('r', r)):
        prove_chain(prime, certificates[name], certificates['leaves'])

    a = 3**h * q
    z = 2**h * v
    rows = [(1, 1, n-2), (a, 2*a, 1), (z, z, 2*z)]
    assert all(min(row) > 0 and sum(row) == n for row in rows)
    expected = [[(1, 514), (1, 1071), (513, 1116)],
                [(1, 513), (513, 167), (1, 1370)],
                [(1, 1), (0, 0), (1, 2)],
                [(1, 1), (1, 1), (0, 0)]]
    depths, common = [], 1
    for p, target in zip((2, 3, q, v), expected):
        actual = [profile(n, row, p) for row in rows]
        assert actual == target
        first, valuations = zip(*actual)
        if all(first):
            depths.append(max(first))
            common *= p**min(valuations)
    assert min(depths) == 513 and common == 2**514 * 3**167
    # Disjointness for d and q, then for each against v, then for r.
    assert d*q == n-1 and v*(2*d+1) > n and v*(2*q+1) > n
    assert 2*r > n and n-r < min(d, q, 2*v)
    for p in (2, 3, q, v):
        value, digit_sum = n, 0
        while value:
            value, digit = divmod(value, p)
            digit_sum += digit
        assert digit_sum >= 3
    print('Weight digits:', len(str(n)), 'sharp depth:', min(depths),
          'minimum number of primes: 4')
    print('Elliptic steps:', *(len(certificates[name]) for name in ('q', 'v', 'r')))


if __name__ == '__main__':
    verify()
