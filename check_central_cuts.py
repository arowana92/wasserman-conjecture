#!/usr/bin/env python3
"""Check the central-cut congruence and enumerate its quadratic equality case."""
import argparse
from collections import Counter
from math import comb, gcd, isqrt, lcm


def exact_congruence(limit):
    pairs = partial = 0
    direct_equality = set()
    for n in range(3, limit + 1):
        coefficients = [comb(n, a) for a in range(n // 2 + 1)]
        for b in range(1, len(coefficients)):
            for a in range(1, b + 1):
                pairs += 1
                common = gcd(coefficients[a], coefficients[b])
                weight_part = gcd(n, common)
                predecessor_part = gcd(n - 1, common)
                h, rem = divmod(weight_part * (n - 2*a) * (n - 2*b), n)
                assert rem == 0
                gap = b - a
                f = (h - weight_part)**2 - 4*(weight_part*gap)**2
                assert f % ((n - 1) // predecessor_part) == 0
                partial += 1 < weight_part < n or 1 < predecessor_part < n - 1
                if weight_part * predecessor_part != 1:
                    continue
                assert gap > 0
                m = f // (n - 1)
                d = n - a - b
                assert m*d*d == h*((h-1)**2 - 4*gap**2 + m) + m*gap**2
                if m == 0:
                    assert h == 2*gap + 1 and 4*b*(n-a) == n*(n-1)
                    direct_equality.add((n, a, b))
    return pairs, partial, direct_equality


def prime_list(limit):
    sieve = bytearray([1]) * (limit + 1)
    sieve[:2] = bytes(2)
    for p in range(2, isqrt(limit) + 1):
        if sieve[p]:
            sieve[p*p::p] = bytes((limit - p*p) // p + 1)
    return [p for p in range(2, limit + 1) if sieve[p]]


def binomial_carry(n, a, p):
    while a:
        if a % p > n % p:
            return True
        a //= p
        n //= p
    return False


def quadratic_candidates(max_gap, max_bits):
    for gap in range(1, max_gap + 1):
        t = 2*gap + 1
        roots = [r for r in range(t) if (r*r - gap*gap) % t == 0]
        for bits in range(2, max_bits + 1):
            power = 1 << bits
            low = max(1, isqrt(t*power) - gap)
            high = isqrt(t*power + gap*gap - 1)
            assert high - low < t
            for residue in roots:
                d = low + (residue - low) % t
                if d > high:
                    continue
                n = (d*d - gap*gap) // t
                a, rem = divmod(n - d - gap, 2)
                b = a + gap
                if rem or not power//2 <= n < power or not 0 < a <= b <= n//2:
                    continue
                if n % 4 != 1 or (a & (n-a) and b & (n-b)):
                    continue
                assert 4*b*(n-a) == n*(n-1)
                yield n, a, b


def factor_candidates(max_gap, max_bits):
    """Independently parametrize the assigned divisors N=u*w at equality."""
    for gap in range(1, max_gap + 1):
        t = 2*gap + 1
        for i in range(1, t + 1, 2):
            j, rem = divmod(t, i)
            if rem or gcd(i, j) != 1:
                continue
            # w*j - u*i = 2*gap, and w is fixed modulo i.
            residue = (-pow(j, -1, i)) % i if i > 1 else 0
            for bits in range(2, max_bits + 1):
                power = 1 << bits
                # Largest w in this progression with u*w < power.
                limit = (gap + isqrt(gap*gap + t*(power-1))) // j
                w = limit - (limit - residue) % i
                u = (w*j - 2*gap) // i
                if u <= j or w <= i or (u*w) % 4 != 1:
                    continue
                n = u*w
                a, b = w*((u-j)//2), u*((w-i)//2)
                assert b-a == gap and 4*b*(n-a) == n*(n-1)
                if not power//2 <= n < power or n-a < power//2:
                    continue
                digits = bin(n).count("1")
                if any(bin(c).count("1") + bin(n-c).count("1") == digits for c in (a, b)):
                    yield n, a, b


def lower_powers_pass(n, a, b):
    return (lcm(gcd(n, a), gcd(n, b)) == n and
            lcm(gcd(n-1, a*(n-a)), gcd(n-1, b*(n-b))) == n-1)


def legendre(n, p):
    valuation = 0
    while n:
        n //= p
        valuation += n
    return valuation


def refinement_masks(n, a, primes):
    """Prime divisors of each possible trinomial refining this cut."""
    masks = set()
    for single in {a, n-a}:
        total = n - single
        valuations = [legendre(n, p) - legendre(single, p) for p in primes]
        for x in range(1, total//2 + 1):
            mask = sum(1 << i for i, p in enumerate(primes)
                       if valuations[i] > legendre(x, p) + legendre(total-x, p))
            masks.add(mask)
    return masks


def endpoint_primes(n):
    primes = []
    for value in (n, n - 1):
        remaining = value
        p = 2
        while p*p <= remaining:
            if remaining % p == 0:
                primes.append(p)
                while remaining % p == 0:
                    remaining //= p
            p += 1 if p == 2 else 2
        if remaining > 1:
            primes.append(remaining)
    return sorted(primes)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gaps", type=int, default=512)
    parser.add_argument("--bits", type=int, default=4096)
    parser.add_argument("--direct", type=int, default=200)
    args = parser.parse_args()
    if not (1 <= args.gaps and 2 <= args.bits and 3 <= args.direct):
        parser.error("positive cut range and valid weight bounds required")
    pairs, partial, exact = exact_congruence(args.direct)
    print(f"exact cut pairs {pairs}, proper intermediate gcd factors {partial}")

    primes = prime_list(997)
    counts = Counter()
    found = set()
    low_power_cases = set()
    for n, a, b in quadratic_candidates(args.gaps, args.bits):
        counts["binary"] += 1
        if not lower_powers_pass(n, a, b):
            continue
        counts["low powers"] += 1
        low_power_cases.add((n, a, b))
        for p in primes:
            if n % p in (0, 1) and binomial_carry(n, a, p) and binomial_carry(n, b, p):
                # Independently verify the witness from factorial valuations.
                vp_n = legendre(n, p)
                assert all(vp_n > legendre(c, p) + legendre(n-c, p) for c in (a, b))
                counts["prime witnesses"] += 1
                break
        else:
            if n > 1_000_000:
                raise RuntimeError(f"unresolved candidate {(n, a, b)}")
            g = gcd(n*(n-1), gcd(comb(n, a), comb(n, b)))
            if g > 1:
                counts["exact witnesses"] += 1
                continue
            found.add((n, a, b))
            ps = endpoint_primes(n)
            first, second = (refinement_masks(n, c, ps) for c in (a, b))
            assert all(u & v for u in first for v in second)
            print(f"binomial gcd one {(n, a, b)}, refinement mask counts {len(first)}, {len(second)}")
    covered_exact = {v for v in exact if v[2]-v[1] <= args.gaps and v[0].bit_length() <= args.bits}
    assert {v for v in found if v[0] <= args.direct} == covered_exact
    independent = {v for v in factor_candidates(args.gaps, args.bits) if lower_powers_pass(*v)}
    assert independent == low_power_cases
    print(dict(counts), "binomial gcd one", len(found))


if __name__ == "__main__":
    main()
