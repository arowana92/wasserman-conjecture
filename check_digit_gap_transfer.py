"""Check the row-scale bound by exact coefficients and transfer by modular arithmetic."""

from collections import Counter
from math import comb, gcd, lcm

from check_weight_transfer import positive_partitions


def primes(n):
    return [i for i in range(2, n + 1) if all(i % j for j in range(2, i))]


def prime_layers(n, r):
    w = r
    while w <= n:
        yield w
        w *= r


def coefficient(row):
    a, b, c = row
    return comb(a + b + c, a) * comb(b + c, b)


def scale_checks():
    checked = Counter()
    for n in range(3, 47):
        rows = [(row, coefficient(row)) for row in positive_partitions(n)]
        for p in (2, 3, 5):
            for q in prime_layers(n - 1, p):
                d, c = divmod(n, q)
                if not c:
                    continue
                for row, k in rows:
                    if k % p == 0:
                        continue
                    high, low = zip(*(divmod(z, q) for z in row))
                    assert sum(high) == d and sum(low) == c
                    if all(d * b == c * x for x, b in zip(high, low)):
                        assert gcd(d, c) >= 3
                        checked["proportional"] += 1
                        continue
                    u = gcd(*row)
                    assert u <= d * c and k % (n // u) == 0
                    for other, h in rows:
                        m = n // gcd(*other)
                        g = gcd(n, k, h)
                        assert u * g % m == 0
                        checked["exact pairs"] += 1
    return checked


def progression(p, d, c, a0, bound):
    n = d * p**a0 + c
    moduli = []
    periods = []
    for r in primes(max(bound * d * c, n)):
        if r == p:
            continue
        w1 = max(prime_layers(bound * d * c, r), default=1)
        w2 = max(prime_layers(n, r), default=1) if n * (n - 1) % r == 0 else 1
        w = w1 * w2
        if w > 1:
            moduli.append(w)
            periods.append((w // r) * (r - 1))
    return lcm(*moduli), lcm(*periods)


def layer(parts, w):
    return sum(z % w for z in parts) // w


def lifted_checks(p, d, c, a0, bound):
    q0 = p**a0
    n = d * q0 + c
    assert 1 <= c < q0
    modulus, period = progression(p, d, c, a0, bound)
    a = a0 + period
    assert pow(p, a, modulus) == q0 % modulus
    tests = Counter()
    for row in positive_partitions(n):
        if coefficient(row) % p == 0:
            continue
        shape = [divmod(z, q0) for z in row]
        if all(d * b == c * x for x, b in shape):
            continue
        for other in positive_partitions(n):
            g0 = gcd(n * (n - 1), coefficient(row), coefficient(other))
            m = n // gcd(*other)
            if m > bound * d * c:
                continue
            direction = [z // (n // m) for z in other]
            for r in primes(n):
                if r == p or n * (n - 1) % r:
                    continue
                v = 0
                for w in prime_layers(n, r):
                    residues = [(x * pow(p, a, w) + b) % w for x, b in shape]
                    # Reconstruct the scaled second row modulo w without inverting m.
                    weight = (d * pow(p, a, m * w) + c) % (m * w)
                    assert weight % m == 0
                    second = [z * (weight // m) for z in direction]
                    assert layer(residues, w) == layer(row, w)
                    assert layer(second, w) == layer(other, w)
                    v += layer(residues, w)
                    tests["matched layer pairs"] += 1
                    if m % r == 0 and n % r == 0:
                        tests["noninvertible layers"] += 1
                valuation = 0
                value = g0
                while value % r == 0:
                    value //= r
                    valuation += 1
                assert v >= valuation
    return tests


if __name__ == "__main__":
    print(scale_checks())
    totals = Counter()
    for args in [(2, 3, 6, 3, 2), (2, 5, 7, 3, 1), (3, 4, 2, 2, 2), (5, 2, 3, 1, 2)]:
        totals.update(lifted_checks(*args))
    print(totals)
