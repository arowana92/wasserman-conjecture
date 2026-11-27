"""Check the finite small-marker cases for 1 + d R p^(hk+j), R=p^h-1."""

from argparse import ArgumentParser
from collections import Counter
from math import comb
from pathlib import Path
from subprocess import run

from flint import fmpz

from check_binary_gap_finite import check_weight, factors
from search_endpoint import splits


MAX_D, MAX_K, MAX_J = 512, 6, 2
LARGE_MARKER = 2**31 - 1


def digit_sum(n, p):
    total = 0
    while n:
        n, digit = divmod(n, p)
        total += digit
    return total


def fits(x, total, p):
    """Each digit of x fits below the corresponding digit of total."""
    if p == 2:
        return x & total == x
    while x:
        x, a = divmod(x, p)
        total, b = divmod(total, p)
        if a > b:
            return False
    return True


def split_rows_possible(n, d, r, p, q, w, counts):
    t = d * r
    for u, _, other_factors in splits(factors(n)):
        v = n // u
        if u <= 2 * r or 2 * u >= t or v <= 2 * w:
            continue
        counts["size-feasible scales"] += 1
        upper = [x for x in range(u, t, u) if fits(x, t, p)]
        if not any(x + y < t and fits(y, t - x, p)
                   for i, x in enumerate(upper) for y in upper[i:]):
            continue
        counts["first row digit splits"] += 1
        bases = [q] + [base for base, _ in other_factors]
        for a in range(1, (u - 1) // (2 * r) + 1):
            for b in range(a, (u - 1) // r - a + 1):
                x, y = r * a * v, r * b * v
                if all(fits(x, n, s) and fits(y, n - x, s) for s in bases):
                    return True
    return False


def prime_power(n):
    if n < 2:
        return None
    powers = factors(n)
    if len(powers) != 1:
        return None
    p = powers[0][0]
    h = 0
    while n > 1:
        n //= p
        h += 1
    return p, h


def markers():
    # Among even R, h>1 requires R=8, p=3, h=2; otherwise j=0.
    # At odd R we have p=2 and j<=2. Thus p^j<=4 throughout.
    bound = (MAX_K + 3) * MAX_D**4 * 4**3 * 2**MAX_K
    candidates = {2**i - delta for i in range(1, bound.bit_length() + 1)
                  for delta in (0, 1)}
    return sorted(r for r in candidates if r <= bound and prime_power(r)
                  and prime_power(r + 1))


def scan_marker(r):
    p, h = prime_power(r + 1)
    q, _ = prime_power(r)
    counts = Counter()
    fallback = set()
    for k in range(1, MAX_K + 1):
        for j in range(min(h - 1, MAX_J) + 1):
            w = p**j * (r + 1)**k
            for d in range(1, MAX_D + 1):
                n = 1 + d * r * w
                counts["parameters"] += 1
                if n <= 10**9 or r > (k + 3) * d**4 * p**(3*j) * 2**k:
                    continue
                if all(s in (p, q) for s, _ in factors(d)):
                    continue
                if digit_sum(d * r, p) > 2 * h * (p - 1):
                    fallback.add(n)
                elif split_rows_possible(n, d, r, p, q, w, counts):
                    fallback.add(n)
    for n in sorted(fallback):
        check_weight(n, counts)
    counts["full CRT weights"] = len(fallback)
    return counts


def elementary_checks():
    for p in (2, 3, 5, 7):
        for n in range(1, 40):
            for x in range(n + 1):
                assert fits(x, n, p) == (comb(n, x) % p != 0)
    r = 2**521 - 1
    assert fmpz(r).is_prime()
    lucas = 4
    for _ in range(519):
        lucas = (lucas * lucas - 2) % r
    assert lucas == 0
    assert r > (48 + 3) * 2**(4*96 + 3*24 + 48)


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--scale-checker", required=True, type=Path,
                        help="compiled check_marker_split_scales.cpp")
    args = parser.parse_args()
    elementary_checks()
    entries = markers()
    assert LARGE_MARKER in entries
    totals = Counter()
    for r in entries:
        if r == LARGE_MARKER:
            result = run([str(args.scale_checker.resolve())], capture_output=True,
                         text=True, check=True)
            assert not result.stdout.strip(), result.stdout
            print(r, result.stderr.strip(), flush=True)
        else:
            counts = scan_marker(r)
            print(r, dict(counts), flush=True)
            totals.update(counts)
    print("Totals excluding the separately enumerated 31-bit marker:", totals)
