"""Certify 1+p^e G(p^a) at every a,e for specified monic quadratic blocks."""

from argparse import ArgumentParser
from collections import Counter
from functools import lru_cache
from itertools import combinations
from math import gcd, isqrt, lcm

from sympy.ntheory.residue_ntheory import sqrt_mod

from check_close_powers import prime_base, subproducts
from check_polynomial_powers import (
    divisor_bound,
    evaluate,
    full_power_bound,
    power_limit,
    row_shapes,
)
from check_six_bit_rays import divisors, normalized_rows
from prime_certificates import free_rows
from search_even_endpoint import avoids, prime_powers


@lru_cache(maxsize=8192)
def root_classes(linear, constant, modulus):
    discriminant = linear**2 - 4 * constant
    roots = sqrt_mod(discriminant, 4 * modulus, all_roots=True)
    return sorted(
        {((z - linear) // 2) % modulus for z in roots if (z - linear) % 2 == 0}
    )


def coprime_part(value, modulus):
    while (shared := gcd(value, modulus)) > 1:
        value //= shared
    return value


def local_obstruction(p, total, parts, assigned, u, companion):
    forced_a = sum(parts[:2]) // gcd(*parts[:2])
    forced_c = lcm(
        sum(companion) // gcd(*companion),
        *(sum(pair) // gcd(*pair) for pair in combinations(companion, 2)),
    )
    candidates = sorted({2, 3, 5, 7} | set(prime_powers(total * u)))
    for q in candidates:
        if q == p:
            continue
        if u * q > 100000:
            continue
        levels = [q]
        while u * levels[-1] * q <= 100000:
            levels.append(levels[-1] * q)
        modulus = u * levels[-1]
        power = 1
        while True:
            n = 1 + power * total
            if n % u == 0:
                v = n // u
                if v % q == 0 and u % q == 0:
                    pass  # U must be unitary in N.
                elif u % q == 0 or assigned % q == 0:
                    actual = [power * x + int(i == 2) for i, x in enumerate(parts)]
                    if forced_a % q and all(
                        sum(x % w for x in actual) < w for w in levels
                    ):
                        break
                elif v % q == 0 or total % q == 0:
                    if forced_c % q and all(
                        sum((v * x) % w for x in companion) < w for w in levels
                    ):
                        break
                else:
                    break
            power = power * p % modulus
            if power == 1:
                return q
    return None


def inspect_row(p, total, parts, assigned, counts):
    x, y, _ = parts
    possible = coprime_part(gcd(x, y), p * total)
    remaining = total // assigned
    if possible < 3 or 27 * remaining**2 > possible**3:
        return
    denominator = total // gcd(*parts)
    if denominator < 364 or any(
        3 * w >= denominator for w in prime_powers(denominator).values()
    ):
        return
    if not avoids(total, parts, p):
        return
    powers = tuple(prime_powers(remaining).values())
    for u in divisors(possible):
        if u < 3 or 27 * remaining**2 > u**3 or any(2 * w + 1 > u for w in powers):
            continue
        # The unknown outer power must satisfy p^e G(X) = -1 modulo U.
        residue = 1
        wanted = -pow(total, -1, u) % u
        while residue != wanted:
            residue = residue * p % u
            if residue == 1:
                break
        if residue != wanted:
            continue
        counts["companion searches"] += 1
        for companion in normalized_rows(u, powers):
            q = local_obstruction(p, total, parts, assigned, u, companion)
            if q is None:
                raise ArithmeticError(
                    ("unresolved assignment", p, total, parts, assigned, u, companion)
                )
            counts["local obstructions"] += 1


def small_inputs(G, counts):
    c, b, _ = G
    for X in range(2, max(b, c) + 1):
        p = prime_base(X)
        if not p:
            continue
        total = evaluate(G, X)
        factors = prime_powers(total)
        for row in free_rows(total, p):
            for shifted in range(3):
                x, y = (v for i, v in enumerate(row) if i != shifted)
                eligible = [w for w in factors.values() if x % w == 0 or y % w == 0]
                counts["overlapping block rows"] += 1
                for assigned in subproducts(eligible):
                    inspect_row(p, total, (x, y, row[shifted]), assigned, counts)


def check_quadratic(b, c, counts):
    G = (c, b, 1)
    low = max(2, b + 1, c + 1)
    small_inputs(G, counts)
    for row in row_shapes(G):
        for shifted in range(3):
            if shifted == 1 and row[0] == row[1]:
                continue
            x, y = (v for i, v in enumerate(row) if i != shifted)
            upper_u = min(
                (x, y, divisor_bound(x, y)), key=lambda f: (len(f), sum(map(abs, f)))
            )
            limits = [full_power_bound(G, divisor_bound(G, f)) for f in (x, y)]
            if any(len(f) != 1 for f in limits):
                raise ValueError(("nonconstant predecessor bound", G, row))
            upper_b = lcm(*(f[0] for f in limits))
            counts["polynomial shapes"] += 1
            for B in divisors(upper_b):
                high = power_limit(G, 1, upper_u, (B,))
                if high is None:
                    raise ArithmeticError(("missing cutoff", G, row))
                if high < low:
                    continue
                factors = tuple(prime_powers(B).values())
                for r in root_classes(b, c, B):
                    if any(evaluate(x, r) % w and evaluate(y, r) % w for w in factors):
                        continue
                    start = low + (r - low) % B
                    for X in range(start, high + 1, B):
                        total = evaluate(G, X)
                        if gcd(total // B, B) != 1:
                            continue
                        values = tuple(evaluate(f, X) for f in (x, y, row[shifted]))
                        if 27 * (total // B) ** 2 > gcd(*values[:2]) ** 3:
                            continue
                        p = prime_base(X)
                        if p:
                            counts["prime-power specializations"] += 1
                            inspect_row(p, total, values, B, counts)


if __name__ == "__main__":
    cli = ArgumentParser(description=__doc__)
    cli.add_argument("--height", type=int, default=10)
    args = cli.parse_args()
    if args.height < 0:
        cli.error("height must be nonnegative")
    counts = Counter()
    for b in range(args.height + 1):
        for c in range(args.height + 1):
            d = b * b - 4 * c
            if d > 0 and isqrt(d) ** 2 == d:
                continue
            check_quadratic(b, c, counts)
            counts["quadratics"] += 1
    print(dict(counts))
