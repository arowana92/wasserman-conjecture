"""Residue and carry certificates for all-exponent prime classes at cofactor 380."""
from fractions import Fraction
from math import gcd, lcm

from carry_blocks import CarryGraph, accepted_hull
from check_cofactor_986_powers import low_prefixes, productive_edges
from prime_certificates import primes_through


SLOPES = (130, 247, 3)
MODULUS = 625 * 169 * 19
PERIOD = 58500


def order(value, modulus):
    residue = 1
    for exponent in range(1, modulus + 1):
        residue = residue * value % modulus
        if residue == 1:
            return exponent
    raise ValueError((value, modulus))


def binary_prefix(w, bits):
    modulus = 4 << bits
    v = (380*w + 1) * pow(13, -1, modulus) % modulus
    return all(v & (v << j) == 0 for j in (1, 2, 3))


def fixed_bases():
    assert accepted_hull(CarryGraph(2, SLOPES, (0, 0, 1))) == (0, 0)
    for p, digit, digits, target in [(5, 3, 4, 458), (13, 4, 2, 56)]:
        edges = productive_edges(p, SLOPES)
        assert {x for x, _ in low_prefixes(edges, p, digit, digits)} == {target}
    edges = productive_edges(19, SLOPES)
    targets = {x for x, _ in low_prefixes(edges, 19, 6, 2)}
    assert targets and pow(6, 19, 361) not in targets
    targets = ((458, 625), (56, 169), (6, 19))
    r = next(x for x in range(458, MODULUS, 625)
             if all(x % m == value for value, m in targets))
    orders = [order(value, m) for value, m in targets]
    assert orders == [500, 78, 9] and lcm(*orders) == PERIOD
    assert lcm(500, 156, 18) == PERIOD
    count = sum(gcd(j, PERIOD) == 1 for j in range(PERIOD))
    assert count == 14400
    # Modulo 16 no admissible W is 1. For the other q == 1 mod 4,
    # changing the exponent by PERIOD preserves h+2 bits of W.
    assert not binary_prefix(1, 4)
    for e in (1, 3):
        classes = []
        for q in range(5, 32, 4):
            h = ((q-1) & (1-q)).bit_length() - 1
            if h < 4 and binary_prefix(q**e, h+2):
                classes.append(q)
        assert len(classes) == 5
        print("Exponent mod 4:", e, "binary base classes mod 32:", classes)
    mass = Fraction(count, 500*156*18) * Fraction(5, 16) / 2 * Fraction(18, 19)
    assert mass == Fraction(3, 1976)
    print("Scale modulo", MODULUS, "=", r, "period", PERIOD)
    return mass


def companion_sieve(limit):
    records = []
    for p in primes_through(limit):
        if (380*PERIOD) % p == 0 or pow(-380, PERIOD, p) != 1:
            continue
        r = -pow(380, -1, p) % p
        for h in range(2, 6):
            m = p**h
            target = pow(r, p**(h-1), m)
            assert target % p == r and pow(target, PERIOD, m) == 1
            v = (380*target + 1) * pow(13, -1, m) % m
            if sum(a*v % m for a in (1, 4, 8)) >= m:
                assert sum((a*target+b) % p for a, b in zip(SLOPES, (0, 0, 1))) >= p
                records.append((p, h))
                break
    return records


if __name__ == "__main__":
    mass = fixed_bases()
    records = companion_sieve(100000)
    print("Companion exclusions (prime, layer):", records)
    for p, h in records:
        mass *= 1 - Fraction(1, (p-1)*p**(h-1))
    assert len(records) == 53 and 1-mass > Fraction(9985366, 10000000)
    print("Certified prime density:", float(1-mass))
