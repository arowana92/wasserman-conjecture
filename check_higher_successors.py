"""Certify N = 1 + p**a * (p**j + d) for j >= 2 and bounded offsets."""

from argparse import ArgumentParser
from collections import Counter
from functools import lru_cache
from itertools import product
from math import gcd, isqrt

from check_carry_intersections import COFACTORS
from check_close_powers import bases
from check_prime_power_cofactors import carry_language
from check_successor_exponents import (
    Case, bounded_exponents, complementary, first_progression, fixed_divisor,
    local_parts, orbit_sieve, powers, required_bases,
)
from prime_certificates import primes_through
from search_endpoint import atoms, avoids, smallest_factors


LAST = Case(2278, 7, 2, (833, 1407, 38), (1, 2, 4), 37, 1369, 5665, 7973)


@lru_cache(maxsize=128)
def logarithms(p, modulus):
    """Map each attainable residue to its least nonnegative exponent."""
    return {v: i for i, v in enumerate(powers(p % modulus, modulus))}


def assignments(Q, p, height, factors):
    d = Q - height
    # The three groups are the powers in e, in f, and in the companion factor T.
    choices = [(0, 1, 2) if r != p else (1,) for r, w in factors]
    for colors in product(*choices):
        e = f = T = 1
        companion_powers = []
        for (_, w), color in zip(factors, colors):
            if color == 0:
                e *= w
            elif color == 1:
                f *= w
            else:
                T *= w
                companion_powers.append(w)
        if e > d:
            continue
        low = 2 * max(companion_powers, default=1) + 1
        upper = (d - 1) // f
        while low <= upper and low**3 < 27*T*T:
            low += 1
        for u in range(low, upper + 1):
            if gcd(u, Q) != 1:
                continue
            step = u * e
            # A shifted small part requires t = x+y divisible by e and t = Q mod u.
            t0 = e * (Q * pow(e, -1, u) % u)
            for y in range(u*f, d, u*f):
                if gcd(e, y) != 1:
                    continue
                for shifted, start in ((False, step), (True, (t0-y) % step or step)):
                    for x in range(start, d-y+1, step):
                        if not shifted and y > x:
                            continue
                        large = Q - x - y
                        first = (y, large, x) if shifted else (y, x, large)
                        if not avoids(first, p):
                            continue
                        progression = first_progression(Q, u, e*f, factors, *first[:2])
                        if progression is None:
                            continue
                        residue, modulus = progression
                        if gcd(p, modulus) != 1:
                            continue
                        orbit = logarithms(p, modulus)
                        if residue in orbit:
                            yield (u, T, tuple(companion_powers), first, residue, modulus,
                                   orbit[residue], len(orbit))


def prime_layer(r, cap):
    w = r
    while w*r <= cap:
        w *= r
    return w


def valid_scale(Q, u, r, w, value):
    power = gcd(u, w)
    return not (1 < power < w and (Q*value + 1) % (r*power) == 0)


def first_prefixes(Q, u, T, first, p, start, period):
    for r in (set(bases(Q//T)) | set(bases(u))) - {p}:
        w = prime_layer(r, 512)
        head = pow(p, start, w)
        step = pow(p, period, w)
        v = head
        x, y, z = first
        while True:
            if valid_scale(Q, u, r, w, v) and avoids((x*v % w, y*v % w, (z*v+1) % w), r):
                break
            v = v*step % w
            if v == head:
                return False
    return True


def local_cases(limit, counts):
    bound = isqrt(limit**5 // 27)
    spf = smallest_factors(bound)
    for p in primes_through(isqrt(bound)):
        height = p*p
        while height < bound:
            for d in range(4, min(limit, bound-height) + 1):
                Q = height + d
                if 27*Q*Q >= d**5:
                    continue
                cofactor = Q
                while cofactor % p == 0:
                    cofactor //= p
                if cofactor <= 1000 and (p == 2 or cofactor not in COFACTORS):
                    continue
                factors = atoms(Q, spf)
                if sum(r != p for r, _ in factors) < 2 or any(w > d for _, w in factors):
                    continue
                for u, T, full, first, R, M, a0, L in assignments(Q, p, height, factors):
                    counts["first rows in the power orbit"] += 1
                    if not first_prefixes(Q, u, T, first, p, a0, L):
                        continue
                    counts["first rows after prefix tests"] += 1
                    for other in complementary(u, full):
                        counts["row pairs"] += 1
                        c = Case(Q, u, T, first, other, p, height, R, M)
                        if orbit_sieve(c, 256) and orbit_sieve(c, 4096):
                            counts["pairs after simultaneous tests"] += 1
                            yield c
            height *= p


@lru_cache(maxsize=8192)
def suffix_possible(r, slopes, state):
    return carry_language(r, slopes, state, range(r))


def extendable(c, r, cap=4096):
    first, _ = required_bases(c)
    is_first = r in first
    slopes = c.first if is_first else c.other
    w = prime_layer(r, cap)
    period = logarithms(c.p, c.modulus)
    start = period[c.residue]
    cycle = powers(c.p % w, w)
    stride = gcd(len(period), len(cycle))
    for i in range(start % stride, len(cycle), stride):
        v = cycle[i]
        if not valid_scale(c.Q, c.u, r, w, v):
            continue
        parts = local_parts(c, v, w, is_first)
        if avoids(tuple(k % w for k in parts), r):
            state = tuple(k // w for k in parts)
            if suffix_possible(r, slopes, state):
                return True
    return False


def solve(limit, exponent_limit, list_local):
    counts = Counter()
    excluded = Counter()
    remaining = set()
    for c in local_cases(limit, counts):
        first, other = required_bases(c)
        r = next((r for r in sorted(set(first + other)) if not extendable(c, r)), None)
        if r:
            excluded[r] += 1
        else:
            r = fixed_divisor(c)
            if r:
                counts["pairs excluded by fixed prime"] += 1
            else:
                remaining.add(c)
        if list_local:
            print(c, "exclusion base", r, flush=True)
    print(dict(counts))
    print("Excluded by nonextendable digits at each base:", dict(sorted(excluded.items())))
    expected = {LAST} if limit >= LAST.Q - LAST.height else set()
    if remaining != expected:
        raise ArithmeticError(("unexpected residual templates", remaining ^ expected))
    if remaining:
        orbit = logarithms(LAST.p, LAST.modulus)
        if (orbit[LAST.residue], len(orbit)) != (28, 48):
            raise ArithmeticError("unexpected exponent class")
        for a in bounded_exponents(LAST, exponent_limit, 7):
            raise ArithmeticError(("unresolved exponent in bounded check", a))
        print(f"Only 1 + 2278 * 37**a, a = 28 mod 48; excluded through a = {exponent_limit}.")
    else:
        print("All exponents excluded from a counterexample.")


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("max_offset", type=int, nargs="?", default=1000)
    parser.add_argument("--exponent-limit", type=int, default=10**24)
    parser.add_argument("--list-local", action="store_true")
    args = parser.parse_args()
    if not 1 <= args.max_offset <= 1000 or args.exponent_limit < 1:
        parser.error("use an offset from 1 to 1000 and a positive exponent limit")
    solve(args.max_offset, args.exponent_limit, args.list_local)
