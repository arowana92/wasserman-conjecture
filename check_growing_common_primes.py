"""Finite threshold, carry and constant checks for growing_common_primes.md."""
from fractions import Fraction as F
from itertools import combinations, combinations_with_replacement
from math import comb, gcd, prod

from check_endpoint_pair_bound import prime_powers, rows
from prime_certificates import free_in_base, primes_through


def constants():
    log_z = F(1, 100)   # in units of sqrt(log X)
    n_tail = 2*F(1, 10)**2/log_z - 2*log_z
    z_tail = 2*F(1, 5)**2/log_z - 2*log_z
    assert n_tail > 1 and z_tail/2 - 2 > 1
    assert F(2, 15) < F(138, 1000) < F(14, 100)
    assert 21*(1-F(1, 25)) > 18+15*F(14, 100)
    for j in (2, 4, 8):
        for e in range(100):
            assert (e+1)**j <= comb(e+2**j-1, 2**j-1)
    print("Digit-tail, weighted-product and divisor-moment comparisons passed")


def two_residue_sieve():
    primes = (3, 5, 7, 11, 13)
    period = prod(primes)
    counts = [0]*(len(primes)+1)
    for n in range(period):
        hits = sum(n % p in (0, 1) for p in primes)
        for j in range(hits+1):
            counts[j] += comb(hits, j)
    remaining = prod(p-2 for p in primes)
    truncated = 0
    for j, actual in enumerate(counts):
        assert actual == sum(2**j*(period//prod(group)) for group in combinations(primes, j))
        truncated += (-1)**j*actual
        assert (truncated-remaining)*(-1)**j >= 0
    assert truncated == remaining
    print(f"All inclusion-exclusion levels passed on {period:,} consecutive residues")


def threshold(masks, primes, k):
    """Minimum over k rows of the largest common listed prime, zero for none."""
    # Equal masks give identical intersections. A strict superset cannot lower
    # the minimum because the subset is available as a replacement row.
    minimal = []
    for mask in sorted(set(masks), key=lambda m: bin(m).count("1")):
        if all(other & mask != other for other in minimal):
            minimal.append(mask)
    common = set(minimal)
    for _ in range(1, k):
        common = {a & b for a in common for b in minimal}
    bit = min(mask.bit_length() for mask in common)
    return primes[bit-1] if bit else 0


def exact_checks(limit=160):
    pair_counts = {p: 0 for p in (3, 5, 11, 17)}
    full_counts = pair_counts.copy()
    total = relation_tests = 0
    for n in range(6, limit+1):
        endpoint = n*(n-1)
        seeds = sorted({p for m in (n, n-1) for p, w in prime_powers(m)})
        r = primes_through(n)[-1]
        candidates = sorted({r} | {p for h in range(6) for p, w in prime_powers(n-h)})
        table = list(rows(n))
        exact_masks = [sum(1 << j for j, p in enumerate(candidates) if value % p == 0)
                       for row, value in table]
        carry_masks = [sum(1 << j for j, p in enumerate(candidates) if not free_in_base(row, p))
                       for row, value in table]
        assert exact_masks == carry_masks
        pair_masks = [sum(1 << j for j, p in enumerate(seeds) if value % p == 0)
                      for row, value in table]
        pair = threshold(pair_masks, seeds, 2)
        full = threshold(exact_masks, candidates, 3)
        assert full <= pair  # the endpoint row is available to the minimizing triple
        for floor in pair_counts:
            pair_counts[floor] += pair > floor
            full_counts[floor] += full > floor

        if n <= 22:
            # Independent calculation with the entire gcd as an integer.
            def largest_factor(number, eligible):
                return max((p for p in eligible if number % p == 0), default=0)
            values = [v for row, v in table]
            assert full == min(largest_factor(gcd(a, b, c), candidates)
                               for a, b, c in combinations_with_replacement(values, 3))
            assert pair == min(largest_factor(gcd(endpoint, a, b), seeds)
                               for a, b in combinations_with_replacement(values, 2))

        for row, value in table:
            a, b, c = row
            assert a*b*value % endpoint == 0
            if n <= 65:
                cut = comb(n, c)
                for other, other_value in table[::max(1, len(table)//11)]:
                    single = gcd(endpoint, other_value)
                    ratio = single//gcd(single, cut)
                    u, s = gcd(n, ratio), gcd(n-1, ratio)
                    assert u*s == ratio
                    for z in (c, n-c):
                        assert z % u == 0 and z*(z-1) % s == 0
                    relation_tests += 1
        total += len(table)
    print(f"{total:,} exact rows for N=6..{limit}, {relation_tests:,} ratio assignments")
    print("Cutoff   weights passing endpoint pairs / full three-row list")
    for floor in pair_counts:
        print(f"{floor:>6}   {pair_counts[floor]:>3} / {full_counts[floor]}")


if __name__ == "__main__":
    constants()
    two_residue_sieve()
    exact_checks()
