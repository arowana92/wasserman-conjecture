"""Check arbitrary and endpoint-only prime certificate minima at two weights."""

import argparse
import subprocess
from itertools import combinations
from math import comb, gcd, isqrt, prod
from pathlib import Path

from flint import fmpz

from check_six_prime_minimum import row_valuation


# n, factors of n, factors of n-1, nearest prime gap, factors of near prime minus 1,
# one endpoint ten-set, number of avoidable fours, number of certifying tens.
CASES = (
    (
        7064355413422079,
        (67, 101, 103, 127, 181, 271, 1627),
        (2, 73, 113, 179, 283, 1129, 7487),
        22,
        (2, 2, 2, 3, 79, 76777, 48529193),
        (103, 113, 127, 179, 181, 271, 283, 1129, 1627, 7487),
        1,
        791,
    ),
    (
        16795978692415487,
        (101, 127, 193, 409, 2447, 6779),
        (2, 53, 71, 97, 691, 5147, 6469),
        30,
        (2, 2, 2, 2, 2, 2411, 4691, 46407983),
        (53, 71, 127, 193, 409, 691, 2447, 5147, 6469, 6779),
        20,
        33,
    ),
)


def require(condition, detail):
    if not condition:
        raise ArithmeticError(detail)


def trial_prime(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))


def prove_prime(n, factors):
    """Lucas primality criterion with the complete factorization of n-1."""
    require(
        prod(factors) == n - 1 and all(trial_prime(p) for p in set(factors)),
        ("prime certificate factors", n),
    )
    for a in range(2, n):
        if pow(a, n - 1, n) == 1 and all(
            gcd(pow(a, (n - 1) // q, n) - 1, n) == 1 for q in set(factors)
        ):
            return a
    raise ArithmeticError(("no generator", n))


def clean_digits(n, a, b, p):
    c = n - a - b
    while n:
        if a % p + b % p + c % p != n % p:
            return False
        n, a, b, c = n // p, a // p, b // p, c // p
    return True


def covered_sixes(n, primes):
    witnessed = set()
    count = 0
    for line in (
        Path(__file__).with_name("eight_prime_witnesses.txt").read_text().splitlines()
    ):
        if not line or line.startswith("#"):
            continue
        weight, a, b = map(int, line.split())
        if weight != n:
            continue
        require(min(a, b, n - a - b) > 0, line)
        avoided = [p for p in primes if row_valuation(n, a, b, p) == 0]
        require(avoided == [p for p in primes if clean_digits(n, a, b, p)], line)
        witnessed.update(combinations(avoided, 3))
        count += 1
    for group in combinations(primes, 6):
        require(
            any(
                t in witnessed and tuple(p for p in group if p not in t) in witnessed
                for t in combinations(group, 3)
            ),
            ("uncovered six", group),
        )
    return count


def four_profile(scanner, n, primes):
    lines = subprocess.check_output(
        [scanner, str(n), *map(str, primes)], text=True
    ).splitlines()
    avoidable = set()
    full = list(combinations(primes, 4))
    require(len(lines) == len(full), "missing scanner output")
    for expected, line in zip(full, lines):
        *group, a, b = map(int, line.split())
        require(tuple(group) == expected, line)
        if (a, b) != (0, 0):
            require(
                0 < a <= b <= n - a - b
                and all(row_valuation(n, a, b, p) == 0 for p in group),
                line,
            )
            avoidable.add(expected)
    return avoidable


def reverse_cuts(n, group):
    """Independently enumerate binomial-unit cuts by decreasing CRT lifts."""
    period = prod(group)
    residues = [0]
    for p in group:
        if n % p:
            m = period // p
            piece = m * pow(m, -1, p)
            residues += [(r + piece) % period for r in residues]
    levels = []
    for p in group:
        power = p * p
        while power <= n:
            levels.append((power, n % power))
            power *= p
    levels.sort(reverse=True)
    cuts = []
    for r in residues:
        a = n // 2 - (n // 2 - r) % period
        while a > 0:
            for power, bound in levels:
                excess = a % power - bound
                if excess > 0:
                    a -= (1 + (excess - 1) // period) * period
                    break
            else:
                cuts.append(a)
                a -= period
    return sorted(cuts)


def independently_avoidable(n, group):
    small = reverse_cuts(n, group)
    cuts = set(small) | {n - a for a in small}
    for i, a in enumerate(small):
        if 3 * a > n:
            break
        for j in range(i, len(small)):
            b = small[j]
            if a + 2 * b > n:
                break
            if a + b in cuts and all(clean_digits(n, a, b, p) for p in group):
                return True
    return False


def two_units(n, p):
    units = []
    power = 1
    while n and len(units) < 2:
        n, digit = divmod(n, p)
        units.extend([power] * min(2 - len(units), digit))
        power *= p
    require(len(units) == 2, ("too few digit units", p))
    return sum(units)


def check(case, scanner, cross_all):
    n, factors, prev, gap, lucas_factors, ten, nfour, nten = case
    primes = tuple(sorted(factors + prev))
    require(
        all(trial_prime(p) for p in primes)
        and len(set(primes)) == len(primes)
        and prod(factors) == n
        and prod(prev) == n - 1,
        "endpoint factors",
    )
    root = prove_prime(n - gap, lucas_factors)
    rows = covered_sixes(n, primes)
    avoidable = four_profile(scanner, n, primes)
    require(len(avoidable) == nfour, avoidable)
    seven = primes[-7:]
    for group in (ten, seven):
        require(not any(q in avoidable for q in combinations(group, 4)), group)
    valid_tens = sum(
        not any(q in avoidable for q in combinations(t, 4))
        for t in combinations(primes, 10)
    )
    require(valid_tens == nten, valid_tens)
    heights = {p: two_units(n, p) for p in primes}
    require(all(heights[p] > gap for p in seven), heights)
    cross_groups = list(combinations(primes if cross_all else seven, 4))
    for q in cross_groups:
        require(
            independently_avoidable(n, q) == (q in avoidable), ("reverse mismatch", q)
        )
    print(
        f"{n}: {rows} supplied rows cover all {comb(len(primes),6)} six-subsets in pairs; "
        f"{nfour}/{comb(len(primes),4)} avoidable four-subsets; {valid_tens} minimum endpoint ten-sets."
    )
    print(
        f"  Reverse enumeration agreed for {len(cross_groups)} subsets. Lucas generator for {n-gap}: {root}."
    )
    print(f"  Minimum unrestricted list: {(*seven,n-gap)}")
    if n == CASES[0][0]:
        require(avoidable == {(2, 179, 1129, 7487)}, avoidable)
        cap = sorted(heights.values(), reverse=True)[6]
        outside_count = list_count = 0
        for g in range(2, cap):
            if not fmpz(n - g).is_prime():
                continue
            eligible = [p for p in primes if heights[p] > g]
            # Every seven chosen here excludes the only jointly avoidable four.
            require(2 not in eligible, (g, eligible))
            outside_count += 1
            list_count += comb(len(eligible), 7)
        require(
            (cap, outside_count, list_count) == (254, 11, 4631),
            (cap, outside_count, list_count),
        )
        print(
            f"  All optimal lists with exterior prime above N/2: {list_count} across {outside_count} primes."
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scanner", type=Path, required=True)
    parser.add_argument(
        "--cross-all", action="store_true", help="reverse-enumerate every four-subset"
    )
    args = parser.parse_args()
    for case in CASES:
        check(case, str(args.scanner.resolve()), args.cross_all)
