"""Check a squarefree endpoint example and the row forced by its prime factors."""

from itertools import permutations
import json
from math import prod
from pathlib import Path

from check_all_prime_depth_witness import profile, prove_leaf
from check_fixed_prime_depth import add_congruence


def crt(residues):
    value, step = 0, 1
    for residue, modulus in residues:
        value, step = add_congruence(value, step, residue, modulus)
    return value, step


def verify():
    h, u, p, big = 2, 13, 3719, 28049188907
    a = 2*p*big
    assert p > u > 10 and big > 2*p+1
    assert (p-1) % u**h == 0 and p % 2**h == 1 + 2**(h-1)
    assert (big+1) % (u*p)**h == 0 and big % 2**(h+1) == 2**h-1
    q0, step = crt([(1, 2), (pow(2, -1, p**h), p**h),
                    (1+big**h, big**(h+1)), ((u-1)*pow(a, -1, u*u), u*u)])
    q = q0 + 168324*step
    n = 1 + a*q
    v, rem = divmod(n, u)
    r = n-6
    assert rem == 0 and v % u == 1 and v > a and q > a
    factors = (2, u, p, big, q, v)
    assert len(set(factors)) == 6 and prod(factors) == n*(n-1)
    certificates = json.loads(Path(__file__).with_name('squarefree_depth_primes.json').read_text())
    for base in (*factors, r):
        prove_leaf(base, certificates)

    first = ((big+1)*q, 2*(p-1)*big*q, (big-1)*q+1)
    other = (v, v, (u-2)*v)
    endpoint = (1, 1, n-2)
    for row in (first, other, endpoint):
        assert sum(row) == n and min(row) > 0
    depth = {}
    for base in factors:
        entries = [profile(n, row, base) for row in (endpoint, first, other)]
        positions = [pos for pos, _ in entries]
        if all(positions):
            depth[base] = max(positions)
    assert min(depth.values()) == h+1 and depth[big] == depth[2] == h+1

    # Solve the units congruences in each position, independently of the formula.
    recovered = set()
    for assignment in permutations((q, p, big)):
        parts = []
        for exception in assignment:
            mods = [(0, u)] + [(int(prime == exception), prime) for prime in assignment]
            value, period = crt(mods)
            assert period > n and 0 < value < n
            parts.append(value)
        assert sum(parts) == n
        recovered.add(tuple(sorted(parts)))
    assert recovered == {tuple(sorted(first))}
    assert profile(n, first, big)[0] == h+1
    assert r > n//2 and r not in factors and n-r < min(u,p,big,q,v)
    print('Weight:', n, 'digits:', len(str(n)))
    print('First common digits by base:', depth)
    print('Universal depth:', h+1, 'with a six-prime list')


if __name__ == '__main__':
    verify()
