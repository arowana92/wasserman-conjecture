"""Check the mixed six-unit argument with direct expansions and exact integers."""

from collections import Counter
from itertools import combinations_with_replacement
from math import comb, gcd, isqrt


def primes(limit):
    return [n for n in range(2, limit + 1)
            if all(n % k for k in range(2, isqrt(n) + 1))]


def units(n, base):
    result = []
    power = 1
    while n:
        n, digit = divmod(n, base)
        result.extend([power] * digit)
        power *= base
    return result


def check_shapes():
    counts = Counter()
    odd = primes(101)[1:]
    for p in odd:
        for d in range(1, 13):
            factor = p**d + 1
            for q in odd:
                if factor % q:
                    continue
                for t in range(1, 14, 2):
                    u = t*d
                    shifted = 1 + p**u
                    first = units(shifted, q)
                    for i in range(max(1, u-2*d), u+2*d+1):
                        second = units(p**i * factor, q)
                        counts['digit comparisons'] += 1
                        if len(first) != 2 or len(second) != 2:
                            continue
                        counts['simultaneous two-unit parts'] += 1
                        v = min(u, i)
                        m = p**(u-v) + 2*p**(i-v)*factor
                        R, B = p**(u-v), p**(i-v)*factor
                        assert B*first[1]-R*second[1] == R*second[0]-B*(first[0]-1)
                        assert (B*first[1]-R*second[1]) % 2
                        assert first[0] <= t*second[0]
                        assert p**u <= t*m*m
                        if t == 1:
                            assert first[0] == first[1]
                            assert (p, q, factor) in {(3, 5, 10), (5, 3, 6)}
                        assert 216*m*m > factor**3
    print(dict(counts))


def check_reciprocal():
    powers = {p**e: (p, e) for p in primes(251)[1:] for e in range(1, 65)}
    links = {(p, q, d, a) for x, (p, d) in powers.items() if (x+1)//2 in powers
             for q, a in [powers[(x+1)//2]]}
    pairs = {(p, q, d, a, f, i) for p, q, d, a in links
             for other, back, f, i in links if (other, back) == (q, p)}
    assert pairs == {(3, 5, 2, 1, 1, 1), (5, 3, 1, 1, 2, 1)}
    print(f'{len(links)} prime-power identities, {len(pairs)} ordered reciprocal pairs.')


def check_coefficients():
    counts = Counter()
    small_primes = primes(200)
    for n in range(3, 201):
        previous = [p for p in small_primes if (n-1) % p == 0 and len(units(n, p)) <= 6]
        weight = [q for q in small_primes if n % q == 0 and len(units(n, q)) <= 6]
        if not previous or not weight:
            continue
        counts['weights'] += 1
        signatures = set()
        for a in range(1, n//3 + 1):
            for b in range(a, (n-a)//2 + 1):
                c = n-a-b
                k = comb(n, a)*comb(n-a, b)
                counts['exact coefficients'] += 1
                signatures.add(gcd(k, n*(n-1)))
        for x, y in combinations_with_replacement(sorted(signatures), 2):
            assert gcd(x, y) > 1
            counts['endpoint gcd pairs'] += 1
    print(dict(counts))


if __name__ == '__main__':
    check_shapes()
    check_reciprocal()
    check_coefficients()
