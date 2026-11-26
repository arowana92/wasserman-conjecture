"""Exact certificate intersections, prescribed gcds and mandatory large primes."""

from functools import reduce
from itertools import combinations, combinations_with_replacement
from math import comb, gcd, prod

from flint import fmpz
from check_scaled_row_blocks import certify_prime
from prime_certificates import primes_through


def coefficients(n):
    for a in range(1, n // 3 + 1):
        first = comb(n, a)
        for b in range(a, (n - a) // 2 + 1):
            yield first * comb(n - a, b)


def check_785():
    n = 785
    primes = primes_through(n)
    core = [p for p in primes if n * (n - 1) % p == 0]
    incompatible = {p: set(primes) - {p} for p in core}
    universal = set(primes)
    for coefficient in coefficients(n):
        universal = {p for p in universal if coefficient % p == 0}
        for p in core:
            if coefficient % p:
                incompatible[p] = {q for q in incompatible[p] if coefficient % q == 0}
    assert not universal
    certificates = set()
    for triple in combinations(core, 3):
        if all(q in incompatible[p] for p, q in combinations(triple, 2)):
            options = sorted(set.intersection(*(incompatible[p] for p in triple)))
            print(triple, options)
            certificates.update(tuple(sorted((*triple, q))) for q in options)
    assert len(certificates) == 33 and not set.intersection(*map(set, certificates))


SMALL = {
    699: ([(2, 41, 233, 349, 691), (3, 139, 233, 349, 691),
           (2, 17, 233, 349, 691), (2, 41, 139, 233, 349)],
          [(1, 1, 697), (3, 48, 648)], [(1, 349, 349), (233, 233, 233)], [233, 349]),
    717: ([(2, 13, 179, 239, 709), (3, 13, 179, 239, 709),
           (2, 5, 179, 239, 709), (2, 13, 179, 239, 701)],
          [(1, 1, 715), (9, 192, 516)], [(239, 239, 239), (1, 179, 537)], [179, 239]),
}


def small_intersections():
    for n, (lists, common, third, targets) in SMALL.items():
        moduli = list(map(prod, lists))
        seen = [set() for _ in moduli]
        for value in coefficients(n):
            for modulus, divisors in zip(moduli, seen):
                divisors.add(gcd(value, modulus))
        assert all(gcd(*group) > 1 for divisors in seen
                   for group in combinations_with_replacement(divisors, 3))
        assert set.intersection(*map(set, lists)) == set(targets)
        for row, target in zip(third, targets):
            rows = [*common, row]
            assert all(min(r) > 0 and sum(r) == n for r in rows)
            assert gcd(*(comb(n, a) * comb(n - a, b) for a, b, _ in rows)) == target
        print(f"Compulsory primes at {n}: {targets}")


def factorial_valuation(n, p):
    total = 0
    while n:
        n //= p
        total += n
    return total


def valuation(row, p):
    assert min(row) > 0
    return factorial_valuation(sum(row), p) - sum(factorial_valuation(x, p) for x in row)


def crt(equations):
    def combine(old, new):
        a, m = old
        b, n = new
        return a + m * ((b - a) * pow(m, -1, n) % n), m * n
    return reduce(combine, equations, (0, 1))


def power_above(p, limit):
    result = 1
    while result <= limit:
        result *= p
    return result


def general_progressions():
    cases = paired = 0
    for d in range(1, 14):
        for h in range(2, 15):
            if gcd(d, h) > 1:
                continue
            powers = [(int(p), int(p**e)) for p, e in fmpz(d * h).factor()]
            a, _ = crt([(1 + w if h % p == 0 else d, p * w) for p, w in powers])
            b = prod(power_above(p, max(a, d * h)) for p, _ in powers)
            assert gcd(a, b) == d
            eqs = []
            for p, w in powers:
                modulus = power_above(p, a + b + d * h)
                factor = gcd(h, modulus)
                modulus //= factor
                target, remainder = divmod(a + b - 1, factor)
                assert remainder == 0
                eqs.append((target * pow(h // factor, -1, modulus) % modulus, modulus))
            q, step = crt(eqs)
            assert gcd(q, step) == 1
            while q <= max(a, b, h) or not fmpz(q).is_prime():
                q += step
            n = h * q + 1
            v, rem = divmod(n, d)
            assert not rem and n > a + b
            if v.bit_length() > 68:
                continue
            primes = [int(s) for s, _ in fmpz(n).factor()]
            primes += [p for p, _ in powers if h % p == 0] + [q]
            endpoint = (1, 1, n - 2)
            first = (a, b, n - a - b)
            second = (1, q, (h - 1) * q)

            def sampled_gcd(third):
                rows = (endpoint, first, third)
                return prod(s**min(valuation(row, s) for row in rows) for s in primes)

            assert sampled_gcd(second) == v
            cases += 1
            if d >= 3 and v > d and fmpz(v).is_prime():
                assert sampled_gcd((v, v, (d - 2) * v)) == q
                paired += 1
    assert cases and paired
    print(f"Checked {cases} coprime progression pairs by factorial floors, {paired} with two prime gcds")


def check_examples():
    exact = 0
    for n in range(447, 50001, 432):
        q = (n - 1) // 2
        if fmpz(q).is_prime():
            row_a = comb(n, 3) * comb(n - 3, 12)
            row_c = n * comb(2 * q, q)
            assert gcd(n * (n - 1), row_a, row_c) == n // 3
            exact += 1
    print(f"Direct binomial products confirmed the gcd formula at {exact} progression weights")

    for p in (102101, 1000000000000000000000003631621):
        n = 3 * p
        q, r = (n - 1) // 2, n - 2
        for prime in (p, q, r):
            certify_prime(prime)
        assert n % 432 == 15 and 2 * 3 * p * q == n * (n - 1)
        e, a, c, d = (1, 1, r), (3, 12, n - 15), (1, q, q), (p, p, p)
        assert all(sum(row) == n for row in (e, a, c, d))
        for triples, target in (((e, a, c), p), ((e, a, d), q), ((e, c, d), 6)):
            actual = prod(s ** min(valuation(row, s) for row in triples) for s in (2, 3, p, q))
            assert actual == target
        assert all(valuation(row, s) == 0 for row, s in ((d, p), (c, q), (e, r)))
        assert all(valuation(row, s) > 0 for row in (e, c, d) for s in (2, 3))
        print(f"N={n}: exact gcds {p}, {q}, 6; fourth prime {r}")


if __name__ == "__main__":
    check_785()
    small_intersections()
    general_progressions()
    check_examples()
