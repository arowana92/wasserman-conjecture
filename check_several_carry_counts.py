"""Exact checks and carry-free search for several_carry_counts.md (Python 3.9+)."""
from fractions import Fraction as F
from itertools import combinations_with_replacement
from math import comb, gcd, prod

from check_binary_windows import ln_bounds
from check_endpoint_pair_bound import prime_powers, rows
from prime_certificates import free_in_base, free_rows


def parameters():
    primes = (2, 3, 5, 7, 11, 13, 17, 19, 23)
    cases = (
        (3, 2, F(1, 16), F(1, 1000), (335, 436), F(9990, 10000), F(1, 10**6)),
        (3, 3, F(1, 64), F(1, 1000), (201, 291, 425, 518, 636, 676, 735, 757, 791),
         F(9988, 10000), F(5, 10**7)),
        (4, 3, F(1, 128), F(1, 10000), (200, 290), F(9997, 10000), F(1, 10**8)),
    )
    for k, t, gamma, eps, numbers, cutoff, eta in cases:
        b = (2-gamma)/k
        for p, number in zip(primes, numbers):
            u = F(number, 1000)
            mu = (p-1)*(F(1, 2)+2*eps)/t
            lo, hi = ln_bounds(F(p))
            entropy = (ln_bounds(sum(u**j for j in range(p)))[1]
                       + mu*ln_bounds(1/u)[1])/lo
            assert entropy + (1-F(1, t))*(1-b) < cutoff < 1-eta
            assert 2*eps**2/hi > eta
    print("All thirteen digit-entropy and tail margins pass rational enclosures")


def digit_sum(n, p):
    total = 0
    while n:
        n, digit = divmod(n, p)
        total += digit
    return total


def exact_assignments():
    checked = partial = 0
    for n in range(4, 21):
        table = list(rows(n))
        table.extend(((a, n-a), comb(n, a)) for a in range(2, n//2+1))
        atoms = [(p, w, h) for h in (0, 1) for p, w in prime_powers(n-h)]
        for row, value in table:
            for p in (2, 3, 5, 7):
                excess = sum(digit_sum(z, p) for z in row)-digit_sum(n, p)
                power = p**(excess//(p-1))
                assert excess % (p-1) == 0 and value % power == 0 and value % (p*power)
            a = max(row)
            assert 2 <= a <= n-2 and value % comb(n, a) == 0
        for k in (3, 4):
            for chosen in combinations_with_replacement(table, k):
                values = [value for row, value in chosen]
                assigned = [1]*k
                has_partial = False
                for p, w, h in atoms:
                    i = min(range(k), key=lambda j: gcd(w, values[j]))
                    retained = w//gcd(w, values[i])
                    if retained > 1:
                        assert sum(z % retained for z in chosen[i][0]) == h
                        has_partial |= retained < w
                    assigned[i] *= retained
                assert prod(assigned)*gcd(n*(n-1), *values) == n*(n-1)
                checked += 1
                partial += has_partial
    cuts = [comb(15, a) for a in (3, 5, 7)]
    assert all(value % 2 for value in cuts)
    assert gcd(*cuts) == 13 and gcd(15*14, *cuts) == 1
    print(f"{checked:,} exact multirow assignments, {partial:,} with partial powers")


def endpoint_masks(n, factors):
    rows_seen = set(free_rows(n, 2)) | set(free_rows(n, 3))
    masks = {sum(1 << j for j, p in enumerate(factors) if not free_in_base(row, p))
             for row in rows_seen}
    if n <= 90:
        exact = {row: value for row, value in rows(n) if value % 6}
        assert set(exact) == rows_seen
        expected = {sum(1 << j for j, p in enumerate(factors) if value % p == 0)
                    for value in exact.values()}
        assert expected == masks
    return masks, len(rows_seen)


def carry_free_range(limit=1500):
    total = 0
    for n in range(3, limit+1):
        factors = [p for m in (n, n-1) for p, w in prime_powers(m)]
        masks, count = endpoint_masks(n, factors)
        total += count
        for a, b, c in combinations_with_replacement(masks, 3):
            assert a & b & c, (n, a, b, c)
    print(f"Through {limit:,}, any three of the {total:,} eligible trinomials "
          "at the same weight share an endpoint prime")


if __name__ == "__main__":
    parameters()
    exact_assignments()
    carry_free_range()
