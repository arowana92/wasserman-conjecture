"""Classify every minimum trinomial prime list at the two weights with minimum eight."""

import argparse
import subprocess
from itertools import combinations
from pathlib import Path

from check_eight_prime_minima import (
    CASES,
    clean_digits,
    four_profile,
    independently_avoidable,
    require,
    row_valuation,
    two_units,
)
from check_heavy_predecessor_crt import factor


def digit_layers(m, p):
    layers = []
    power = 1
    while m:
        m, d = divmod(m, p)
        layers.append((power, d))
        power *= p
    return tuple(reversed(layers))


def successor(x, p, layers, end):
    """Next number >= x whose digits fit; x and end are at most the digit target."""
    higher = end + 1
    for power, allowed in layers:
        place = x // power
        digit = place % p
        if digit > allowed:
            return higher
        if digit < allowed:
            higher = (place + 1) * power
    return x


def shared_cuts(m, p, r, low, high, fuel):
    left, right = digit_layers(m, p), digit_layers(m, r)
    x = low
    while x <= high:
        if fuel is not None:
            if fuel[0] == 0:
                return
            fuel[0] -= 1
        a = successor(x, r, right, high)
        if a > high:
            break
        x = successor(a, p, left, high)
        if x == a:
            yield x
            x += 1


def digit_witness(n, p, r, limit=None):
    # Only a returned, independently checked witness is used to eliminate a prime.
    fuel = None if limit is None else [limit]
    for a in shared_cuts(n, p, r, 1, n // 3, fuel):
        for b in shared_cuts(n - a, p, r, a, (n - a) // 2, fuel):
            return a, b
    return None


def two_digit_witness(n, p, r):
    """Exhaust all common small parts when N has two base-r digits."""
    h, j = divmod(n, r)
    require(h < r, ("requires two digits", n, r))
    cuts = []
    for i in range(h // 2 + 1):
        block = i * r
        for delta in range(n % p + 1):
            for y in range((delta - block) % p, j + 1, p):
                a = block + y
                if 0 < a <= n // 2 and clean_digits(n, a, 0, p):
                    cuts.append(a)
    cuts.sort()
    full = set(cuts) | {n - a for a in cuts}
    for i, a in enumerate(cuts):
        if 3 * a > n:
            break
        for k in range(i, len(cuts)):
            b = cuts[k]
            if a + 2 * b > n:
                break
            if a + b in full and a % r + b % r <= j and clean_digits(n, a, b, p):
                return a, b
    return None


def endpoint_options(n, primes):
    heights = sorted((two_units(n, p) for p in primes), reverse=True)
    candidates = sorted(
        {r for j in range(heights[6]) for r, _ in factor(n - j)} - set(primes)
    )
    options = {}
    for r in candidates:
        eligible = set(primes)
        for p in primes:
            witness = digit_witness(n, p, r, 30000)
            if witness:
                require(
                    all(row_valuation(n, *witness, q) == 0 for q in (p, r)),
                    (n, r, p, witness),
                )
                eligible.remove(p)
                if len(eligible) < 7:
                    break
        if len(eligible) < 7:
            continue
        # The surviving primes in these data all exceed N/7850. Exhaustion here
        # certifies all negative pairs, independently of the bounded search.
        for p in sorted(eligible):
            witness = two_digit_witness(n, p, r)
            if witness:
                require(
                    all(row_valuation(n, *witness, q) == 0 for q in (p, r)),
                    (n, r, p, witness),
                )
                eligible.remove(p)
            else:
                # Recheck each negative by an unbounded traversal of digit ideals.
                require(
                    not digit_witness(n, p, r),
                    ("digit traversal disagreement", n, p, r),
                )
        if len(eligible) >= 7:
            options[r] = tuple(sorted(eligible))
    return len(candidates), options


def scan_groups(scanner, n, groups):
    groups = sorted(groups)
    request = "".join(" ".join(map(str, (n, *q))) + "\n" for q in groups)
    output = subprocess.check_output(
        [scanner, "--groups"], input=request, text=True
    ).splitlines()
    require(len(output) == len(groups), "missing groups")
    avoided = set()
    for group, line in zip(groups, output):
        *q, a, b = map(int, line.split())
        require(tuple(q) == group, line)
        if a or b:
            require(
                0 < a <= b <= n - a - b
                and all(row_valuation(n, a, b, p) == 0 for p in q),
                line,
            )
            avoided.add(group)
        else:
            require(
                not independently_avoidable(n, group),
                ("reverse group traversal", n, group),
            )
    return avoided


def classify(n, primes, scanner, dest):
    candidates, options = endpoint_options(n, primes)
    sevens = {q for eligible in options.values() for q in combinations(eligible, 7)}
    fours = four_profile(scanner, n, primes)
    fives = {
        f
        for s in sevens
        for f in combinations(s, 5)
        if all(q in fours for q in combinations(f, 4))
    }
    require(
        not scan_groups(scanner, n, fives),
        ("row avoids five primes of a candidate list", n),
    )
    # With no avoidable five, a pair can avoid seven only by leaving a
    # complementary triple to the companion of a row avoiding four.
    complements = {
        s: [tuple(p for p in s if p not in q) for q in combinations(s, 4) if q in fours]
        for s in sevens
    }
    triples = scan_groups(scanner, n, {t for cs in complements.values() for t in cs})
    good = {s for s, cs in complements.items() if triples.isdisjoint(cs)}
    counts = {}
    for r, eligible in options.items():
        valid = [s for s in combinations(eligible, 7) if s in good]
        if valid:
            counts[r] = len(valid)
            if dest:
                for s in valid:
                    print(n, *s, r, file=dest)
    near = [v for r, v in counts.items() if 2 * r > n]
    actual = (
        candidates,
        len(counts),
        sum(counts.values()),
        min(counts),
        len(near),
        sum(near),
    )
    expected = (
        (520, 50, 24589, 900032540887, 11, 4631)
        if n == CASES[0][0]
        else (799, 59, 3463, 2382070442833, 13, 634)
    )
    require(actual == expected, (n, actual, expected))
    print(
        f"{n}: {candidates} exterior candidates, {len(counts)} occur, "
        f"{sum(counts.values())} minimum lists. Smallest possible largest prime {min(counts)}. "
        f"Above N/2: {len(near)} exterior primes and {sum(near)} lists.",
        flush=True,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scanner", required=True, type=Path)
    parser.add_argument(
        "--lists", type=Path, help="write all lists, one per line after the weight"
    )
    args = parser.parse_args()
    output = args.lists.open("w") if args.lists else None
    try:
        for n, fs, pre, *_ in CASES:
            classify(n, tuple(sorted(fs + pre)), str(args.scanner.resolve()), output)
    finally:
        if output:
            output.close()
