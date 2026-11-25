"""Search for an even-weight endpoint counterexample whose other row covers all odd primes."""

from argparse import ArgumentParser
from collections import Counter
from math import isqrt, prod
from flint import fmpz


def prime_sieve(limit):
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[:2] = b"\x00\x00"
    for p in range(2, isqrt(limit) + 1):
        if sieve[p]:
            count = (limit - p*p) // p + 1
            sieve[p*p::p] = b"\x00" * count
    return [p for p in range(3, limit + 1, 2) if sieve[p]]


def prime_powers(n):
    factors = {int(p): int(p)**int(e) for p, e in fmpz(n).factor()}
    if prod(factors.values()) != n or not all(fmpz(p).is_prime() for p in factors):
        raise ArithmeticError(("factorization", n))
    return factors


def avoids(n, row, p):
    w = p
    while w <= n:
        if n // w != sum(x // w for x in row):
            return False
        w *= p
    return True


def candidate_rows(n, m, factors):
    # A subset of prime powers specifies all the exceptional positions for one part.
    modulus = n - 1
    residues = [0]
    for w in factors.values():
        unit = (modulus//w) * pow(modulus//w, -1, w)
        residues.extend([(r + unit) % modulus for r in residues])
    full_mask = len(residues) - 1
    good = {mask for mask, part in enumerate(residues)
            if 0 < mask < full_mask and part % m == 0}
    for first in good:
        for second in good:
            if second <= first or first & second:
                continue
            third = full_mask ^ first ^ second
            if third <= second or third not in good:
                continue
            row = residues[first], residues[second], residues[third]
            if sum(row) == n:
                yield row


def search(bits, sieve_limit=3000000, verbose=False):
    if not 3 <= bits <= 64 or sieve_limit < 3:
        raise ValueError("use a bit bound from 3 to 64 and a sieve limit of at least 3")
    primes = prime_sieve(sieve_limit)
    counts = Counter()
    for k in range(2, bits):
        power = 1 << k
        last = min(((1 << bits) - 1) // power, isqrt(power // 27) + 1)
        if last < 1:
            continue
        # Entry j represents the odd multiplier m = 2*j + 1.
        remaining = [power*(2*j + 1) - 1 for j in range((last + 1)//2)]
        for q in primes:
            start = ((pow(power, -1, q) - 1) * ((q + 1)//2)) % q
            for index in range(start, len(remaining), q):
                value = remaining[index] // q
                while value % q == 0:
                    value //= q
                remaining[index] = value
        for index, residual in enumerate(remaining):
            m = 2*index + 1
            n = power*m
            modulus = n - 1
            if 27*modulus**2 > power**3:
                continue
            counts["weights"] += 1
            ceiling = modulus // (m + 1)**2
            known = modulus // residual
            if residual > 1:
                # Every prime of the residual exceeds the sieve limit. If too
                # few groups can exceed m, the three exceptional bins cannot fit.
                if residual <= sieve_limit**2:
                    if residual > ceiling or known <= m*m:
                        continue
                    residual_factors = {residual: residual}
                else:
                    if residual <= sieve_limit**3 and (known <= m or residual > ceiling**2):
                        continue
                    residual_factors = prime_powers(residual)
            else:
                residual_factors = {}
            factors = {**prime_powers(known), **residual_factors}
            if len(factors) < 3 or max(factors.values()) > ceiling:
                continue
            for row in candidate_rows(n, m, factors):
                counts["crt_rows"] += 1
                if bin(m).count("1") >= 3:
                    counts["crt_at_binary_eligible_weights"] += 1
                if not all(avoids(n, row, p) for p in factors):
                    continue
                counts["shift_free"] += 1
                if bin(m).count("1") >= 3 and all(avoids(n, row, p) for p in prime_powers(m)):
                    counts["counterexamples"] += 1
                    print("Candidate requiring independent verification", n, row, flush=True)
        if verbose:
            print(f"through k={k}: {dict(counts)}", flush=True)
    return counts


if __name__ == "__main__":
    cli = ArgumentParser(description=__doc__)
    cli.add_argument("bits", type=int)
    cli.add_argument("--verbose", action="store_true")
    args = cli.parse_args()
    print(dict(search(args.bits, verbose=args.verbose)))
