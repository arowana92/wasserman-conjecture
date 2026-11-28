"""Verify a 728-digit weight, the carry profiles, and its all-rows conclusion."""

import gzip
import json
from itertools import permutations
from math import gcd, isqrt, prod
from pathlib import Path

from check_fixed_prime_depth import D, SCALARS, add_congruence, first_row


def curve_sum(p, a, first, second):
    if first is None:
        return second
    if second is None:
        return first
    x, y = first
    u, v = second
    if x == u:
        if (y + v) % p == 0:
            return None
        assert y == v
        slope = (3*x*x + a) * pow(2*y, -1, p) % p
    else:
        slope = (v-y) * pow(u-x, -1, p) % p
    z = (slope*slope-x-u) % p
    return z, (slope*(x-z)-y) % p


def curve_multiple(n, a, point, m):
    acc = None
    while m:
        if m % 2:
            acc = curve_sum(n, a, acc, point)
        m //= 2
        if m:
            point = curve_sum(n, a, point, point)
    return acc


def prove_leaf(n, small):
    if n < 10**6:
        assert n >= 2 and all(n % p for p in range(2, isqrt(n)+1))
        return
    base, factors = small[str(n)]
    assert prod(p**e for p, e in factors) == n-1
    assert pow(base, n-1, n) == 1
    for p, e in factors:
        assert 1 < p < n and e > 0
        prove_leaf(p, small)
        assert gcd(n, pow(base, (n-1)//p, n)-1) == 1


def prove_chain(n, chain, small):
    for current, trace, scale, a, point in chain:
        assert current == n and n > 3 and gcd(n, 6) == 1
        m = n + 1 - trace
        q, rem = divmod(m, scale)
        assert scale > 0 and not rem and trace*trace < 4*n and 1 < q < n
        # q > (n^(1/4)+1)^2, equivalently (sqrt(q)-1)^4 > n.
        left = q*q + 6*q + 1 - n
        assert left > 0 and left**2 > 16*q*(q+1)**2
        x, y = (t % n for t in point)
        b = (y*y - x**3 - a*x) % n
        assert gcd(4*a**3 + 27*b*b, n) == 1
        multiple = curve_multiple(n, a, (x,y), scale)
        assert multiple is not None
        assert curve_multiple(n, a, multiple, q) is None
        n = q
    prove_leaf(n, small)


def profile(n, row, p):
    q = p
    layer = 0
    first = total = carry = 0
    parts = row
    while q <= n:
        layer += 1
        carry = (carry + sum(x % p for x in parts)) // p
        parts = tuple(x//p for x in parts)
        assert carry == n//q - sum(parts)
        if carry and first == 0:
            first = layer
        total += carry
        q *= p
    return first, total


def suffix_layers(p, h, seed, suffix, extra):
    # Local value of q=W/2 and of the weight in the shifted prefix.
    w = seed + p**h * suffix
    mod = p**(h + extra)
    q = w * pow(2, -1, mod) % mod
    n = (1 + D*w) % mod
    cuts = list(range(1, 52154))
    counts = []
    for j in range(1, h + extra + 1):
        m = p**j
        cuts = [k for k in cuts if k*(q % m) % m <= n % m]
        if j >= h:
            counts.append(len(cuts))
    return cuts, counts


def verify():
    w0, modulus = 0, 1
    for p, w in SCALARS.items():
        w0, modulus = add_congruence(w0, modulus, w, p**129)
    assert modulus == 365078**129
    w = w0 + 704588 * modulus
    q = w // 2
    n = 1 + D*w
    v = n // 7
    assert n == 7*v and n-1 == 2*89*293*q and len(str(n)) == 728
    with gzip.open(Path(__file__).with_name('all_prime_depth_certificates.json.gz'), 'rt') as stream:
        certificates = json.load(stream)
    big, rem = divmod(n-24, 151)
    assert rem == 0 and big > 151
    for p in (2,7,89,293):
        assert all(p % i for i in range(2, isqrt(p)+1))
    for p, key in ((v, 'v'), (q, 'q'), (big, 'r')):
        prove_chain(p, certificates[key], certificates['leaves'])
    factors = [2, 7, 89, 293, q, v]
    assert len(set(factors)) == 6 and prod(factors) == n*(n-1)
    a, c = first_row(w), (v, 3*v, 3*v)
    assert sum(a) == sum(c) == n and min(*a, *c) > 0
    expected = [(138,2267,1,2584), (131,734,1,739), (130,237,1,430),
                (130,172,1,339), (0,0,1,2), (1,1,0,0)]
    common = []
    for p, entry in zip(factors, expected):
        first, other = profile(n, a, p), profile(n, c, p)
        assert (*first, *other) == entry
        if min(first[1], other[1]) > 0:
            common.append((p, max(first[0], other[0])))
    assert prod(p for p, _ in common) == 365078
    assert min(j for _, j in common) == 130

    # Enumerate the row forced to take 7 and each predecessor factor except 2.
    powers = [89,293,q]
    shift = prod(powers)
    assert shift > v and q > 7
    rows = set()
    for order in permutations(range(3)):
        x = [0,0,0]
        for pos, p in zip(order, powers):
            other = shift//p
            x[pos] = other * ((v*pow(other, -1, p)) % p)
        if min(x) > 0 and sum(x) == v:
            rows.add(tuple(sorted(7*t for t in x)))
    assert rows == {tuple(sorted(a))}
    # The four avoidance sets indexed by v, q, 293, big are pairwise disjoint.
    assert v > 7 and q > 7 and n % q == n % 293 == 1
    possible = list(range(1, 52154))
    layer = 0
    modulus = 1
    while possible:
        modulus *= 293
        layer += 1
        assert modulus <= n
        possible = [k for k in possible if k*(q % modulus) % modulus <= n % modulus]
    assert layer == 137
    example = (9669*q, 42485*q, 1)
    assert sum(example) == n
    assert profile(n, example, q)[1] == 0 and profile(n, example, 293)[0] == 137

    # Adding 89 reduces the depth sufficient for all triples to the sharp 130.
    possible = list(range(1, 52154))
    for p in (89, 293):
        modulus = 1
        for j in range(130):
            modulus *= p
            possible = [k for k in possible if k*(q % modulus) % modulus <= n % modulus]
    assert not possible
    for modulus, base in ((v,v), (q,q), (293**129,293), (89**129,89)):
        target = n % modulus
        assert gcd(151, target-24) == 1
        assert 2*151*(target+24) < modulus
        assert base > 24 and n % base <= 1

    for p, seed, suffix, length in ((89,9587,75,1), (293,1654,2726850956040050170,8)):
        assert ((w - seed) // p**129) % p**length == suffix
    for h in (16, 129):
        a89, counts = suffix_layers(89, h, 9587, 75, 1)
        assert counts == [254,128]
        a293, _ = suffix_layers(293, h, 1654, 2726850956040050170, 1)
        assert set(a89).isdisjoint(a293)
        _, counts = suffix_layers(293, h, 1654, 2726850956040050170, 8)
        assert counts == [120,60,16,14,8,8,2,2,0]

    # A universal prime would divide the endpoint and have digit sum at most two.
    for p in factors:
        count = 0
        value = n
        while value:
            value, digit = divmod(value, p)
            count += digit
        assert count >= 3
    print('Weight digits:', len(str(n)), 'first common layer: 130, endpoint gcd: 365078')
    print('Minimum prime list size: 4. Last q-cut removed at 293-layer', layer)
    print('Elliptic certificates:', *(len(certificates[x]) for x in ('q','v','r')))


if __name__ == '__main__':
    verify()
