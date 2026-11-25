"""Check the remaining odd cofactor shapes using the binary digits of V."""

from argparse import ArgumentParser
from collections import defaultdict
from functools import lru_cache
from math import gcd, lcm

from search_endpoint import smallest_factors


# d, the three multipliers of Q in A, and the two prime factors of d.
CASES = ((589, (133, 434, 22), (19, 31)),
         (629, (259, 357, 13), (17, 37)),
         (703, (133, 518, 52), (19, 37)))


@lru_cache(None)
def separated(length):
    """All length-bit numbers with nonzero bits at least 3 places apart."""
    if length <= 0:
        return (0,)
    top = 1 << (length - 1)
    return separated(length - 1) + tuple(top + x for x in separated(length - 3))


def candidates(d, a, b, bits):
    # Q=(7V-1)/d is fixed modulo d because A must avoid both primes of d.
    modulus = d*d
    residue = (1 + d*pow(a+b, -1, d))*pow(7, -1, modulus) % modulus
    split = bits // 2
    lows = defaultdict(list)
    for low in separated(split):
        lows[low % modulus].append(low)
    for high in separated(bits - split):
        start = high << split
        for low in lows.get((residue - start) % modulus, ()):
            v = start + low
            if v & ((v << 1) | (v << 2)) or v % 7 == 0:
                continue
            q = (7*v - 1) // d
            if q > 1 and q % 2 and gcd(d, q) == 1:
                yield v, q


def factorial_valuation(n, p):
    value = 0
    while n:
        n //= p
        value += n
    return value


def check_case(d, row, primes, bits):
    a, b, c = row
    assert a+b+c == d and primes[0]*primes[1] == d
    tested = 0
    for v, q in candidates(d, a, b, bits):
        n = d*q+1
        assert n == 7*v and (a+b)*q % d == 1
        actual = (a*q, b*q, c*q+1)
        for p in (7,) + primes:
            if factorial_valuation(n, p) > sum(factorial_valuation(x, p) for x in actual):
                break
        else:
            raise AssertionError(f"Needs further checking: d={d}, Q={q}, V={v}")
        tested += 1
    return tested


def check_generation():
    # Independent direct integer enumeration and digit addition at short lengths.
    for bits in range(12):
        brute = {v for v in range(1 << bits)
                 if all(((v >> j) & 7) in (0, 1, 2, 4) for j in range(bits))}
        assert set(separated(bits)) == brute
    for d, (a, b, _), _ in CASES:
        brute = set()
        for v in separated(26):
            q, rem = divmod(7*v - 1, d)
            if not rem and q > 1 and q % 2 and v % 7 and gcd(d, q) == 1 and (a+b)*q % d == 1:
                brute.add((v, q))
        assert set(candidates(d, a, b, 26)) == brute


def check_base_classes():
    """Check the unit groups and the additional exclusion for small prime bases."""
    spf = smallest_factors(2896)
    small = []
    for d, (a, b, c), (p1, p2) in CASES:
        modulus = 7*d
        target = next(r for r in range(modulus) if (d*r+1) % 7 == 0 and (a+b)*r % d == 1)
        exponent = lcm(6, p1-1, p2-1)
        order = next(k for k in range(1, exponent+1) if pow(target, k, modulus) == 1)
        units = [v for v in range(1, modulus) if gcd(v, modulus) == 1]
        classes = {pow(target, k, modulus) for k in range(order) if gcd(k, order) == 1}
        assert all(pow(v, 2, modulus) != target and pow(v, 3, modulus) != target for v in units)
        direct = {v for v in units if any(pow(v, k, modulus) == target for k in range(1, exponent+1))}
        assert direct == classes
        for q in range(2, len(spf)):
            if spf[q] != q or q % modulus not in classes:
                continue
            if factorial_valuation(d, q) > sum(factorial_valuation(k, q) for k in (a, b, c)):
                continue
            small.append((d, q))
            # Every possible exponent is odd. In four binary digits an
            # interior run of one or two ones cannot belong to a 7V row.
            residues = {(d*pow(q, k, 16)+1) % 16 for k in (1, 3)}
            possible = {7*v % 16 for v in separated(4)}
            assert residues.isdisjoint(possible)
        print(f"d={d}: {len(classes)}/{len(units)} unit classes, target order {order}")
    assert small == [(589, 937), (589, 1973), (629, 2129)]


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("bits", nargs="?", type=int, default=64, choices=range(1, 81),
                        metavar="BITS")
    bits = parser.parse_args().bits
    check_generation()
    check_base_classes()
    for d, row, primes in CASES:
        print(f"d={d}: {check_case(d, row, primes, bits)} candidates excluded", flush=True)
