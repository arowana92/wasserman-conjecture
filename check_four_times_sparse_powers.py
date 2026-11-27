"""Check the exceptional binary cuts at a=4 and exclude composite powers by base."""

from argparse import ArgumentParser
from math import comb

from prime_certificates import primes_through


def odd_cut(b):
    return b % 3 == 2 and ((b + 1) & (3 * b)) == 0


def generated(bits):
    for j in range(bits // 6 + 1):
        odd_part = ((1 << (6 * j + 3)) + 1) // 3
        for h in range(2, bits - 6 * j - 1):
            yield (odd_part << h) - 1


def valuation_two(n):
    return (n & -n).bit_length() - 1


def log_at_two(base, target, h, precision):
    """Recover the exponent modulo 2^precision in the group 1 mod 2^h."""
    modulus = 1 << (h + precision)
    residue, product, step = 0, 1, base % modulus
    for i in range(precision):
        if (product - target) % (1 << (h + i + 1)):
            product = product * step % modulus
            residue += 1 << i
        step = step * step % modulus
    # Independent modular exponentiation checks the answer and its uniqueness.
    assert pow(base, residue, modulus) == target % modulus
    assert pow(base, 1 << (precision - 1), modulus) != 1
    assert pow(base, 1 << precision, modulus) == 1
    return residue


def check(limit, bits):
    for nbits in range(1, bits + 1):
        candidates = list(generated(nbits))
        assert len(set(candidates)) == len(candidates)
        assert len(candidates) == sum(max(0, nbits - 6 * j - 3)
                                      for j in range(nbits // 6 + 1))
    assert sorted(generated(bits)) == [b for b in range(1, 1 << bits, 2) if odd_cut(b)]
    for b in range(1, 300, 2):
        assert odd_cut(b) == (b % 3 == 2 and comb(4*b + 1, b + 1) % 2 == 1)
    for q in (11, 23, 47):
        h = valuation_two(q + 1)
        for e in (0, 1, 3, 128, 65537, 10**16 - 1):
            target = pow(-q, e, 1 << (h + 128))
            assert log_at_two(-q, target, h, 128) == e

    checked = 0
    min_large = None
    survivors = []
    precision = 128
    bound = 10**16
    for q in primes_through(limit):
        if q % 12 != 11:
            continue
        checked += 1
        h = valuation_two(q + 1)
        for m in range(3, precision, 2):
            if odd_cut(q**m):
                survivors.append((q, m))
        modulus = 1 << (h + precision)
        target = 1 - (1 << h) * pow(3, -1, modulus)
        exponent = log_at_two(-q, target, h, precision)
        min_large = exponent if min_large is None else min(min_large, exponent)
        if precision <= exponent < bound and exponent % 2:
            # Keep a possible solution for further checking without constructing q^m.
            survivors.append((q, exponent))
    print(f"binary cuts below 2^{bits}: {len(list(generated(bits)))}")
    print(f"eligible prime bases through {limit}: {checked}")
    print(f"minimum lifted exponent residue: {min_large}")
    print(f"unexcluded composite powers: {survivors}")
    assert not survivors


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--limit", type=int, default=1_000_000)
    parser.add_argument("--bits", type=int, default=20)
    args = parser.parse_args()
    if not (args.bits >= 1 and args.limit >= 2):
        parser.error("need a positive digit range and a prime limit of at least two")
    check(args.limit, args.bits)
