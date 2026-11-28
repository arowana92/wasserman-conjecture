"""Verify ternary factor witnesses with factorial floors and CRT idempotents."""

from argparse import ArgumentParser
from collections import Counter
from itertools import combinations, combinations_with_replacement
from math import gcd, isqrt, lcm, prod
from pathlib import Path

from flint import fmpz


def prime_power(number, prime):
    if not fmpz(prime).is_prime() or number % prime:
        raise ValueError(("invalid prime", number, prime))
    rest = number
    while rest % prime == 0:
        rest //= prime
    return number // rest


def carry(a, b, modulus):
    return (a + b) // modulus > a // modulus + b // modulus


def verify(index, factors, stats):
    X = 3**index
    s, D, M = X + 1, X + 2, X + 4
    # n=3 represents the two relevant A parts modulo q^(2v) for q^v || s
    # and the two C parts vanishing modulo each prime power of M.
    Q = X**3 // 3
    scale = (1 + Q * D * M) // s
    first, second = 3 * Q * s, Q * s * D
    u_cap, t_min = s, M
    for q in factors[0]:
        w = prime_power(s, q)
        level = w * q
        while level <= w * w and not carry(first, second, level):
            level *= q
        if level > w * w:
            raise ArithmeticError(("weight prime", index, q))
        u_cap //= w
    for q in factors[1]:
        w = prime_power(D, q)
        odd = q - 1
        while odd % 2 == 0:
            odd //= 2
        if pow(2, odd, q) == 1:
            raise ArithmeticError(("predecessor prime", index, q))
        t_min *= w
    if u_cap**3 < 27 * t_min**2:
        stats["size"] += 1
        return
    ws = [prime_power(M, q) for q in factors[2]]
    if prod(ws) != M:
        raise ArithmeticError(("incomplete modulus", index))
    # Build each allowed coordinate as a sum of orthogonal idempotents.
    residues = {0}
    for w in ws:
        m = M // w
        increment = s * m * pow(m, -1, w)
        residues |= {(z + increment) % M for z in residues}
    allowed = {z for z in residues if 0 < z < s and z % (s // u_cap) == 0}
    for x, y in combinations_with_replacement(sorted(allowed), 2):
        z = s - x - y
        if z < y or z not in allowed:
            continue
        row = (x, y, z)
        if any(sum(v % w != 0 for v in row) != 1 for w in ws):
            continue
        stats["candidate rows"] += 1
        if any(gcd(M, (a + b) // gcd(a, b)) > 1 for a, b in combinations(row, 2)):
            continue
        for w in ws:
            a, b = [v * scale for v in row if v % w == 0]
            if carry(a, b, w * w):
                break
        else:
            raise ArithmeticError(("unresolved row", index, row))
    stats["CRT"] += 1


def verify_outer_limit(limit, last):
    # n > a/2 since 2^(2n) > X=3^a. Work down in a so that all
    # consecutive LCMs can be built incrementally without factorization.
    multiple = 1
    highest = 0
    count = 0
    for a in range(isqrt(2 * (limit + 1)), last, -1):
        X = 3**a
        F, M, s = X + 2, X + 4, X + 1
        max_n = (limit + 1) // a
        if max_n % 2 == 0:
            max_n -= 1
        if max_n < 1:
            continue
        upper = s // (s & -s)
        threshold = 27 * (F * M) ** 2
        if upper**3 * ((1 << max_n) - 1) ** 2 < threshold:
            continue
        for n in range(highest + 1, max_n + 1):
            multiple = lcm(multiple, n)
        highest = max_n
        # Initially allow even orders too. Remove the 2-part only if needed.
        for exp in (multiple, multiple // (multiple & -multiple)):
            cap = gcd(F, int(pow(fmpz(2), fmpz(exp), fmpz(F))) - 1)
            if upper**3 * cap**2 < threshold:
                break
        else:
            raise ArithmeticError(("outer exponent bound", a, max_n))
        count += 1
    print(f"Verified every outer exponent through {limit} with {count} modulus checks.")


if __name__ == "__main__":
    cli = ArgumentParser(description=__doc__)
    cli.add_argument("--outer-limit", type=int)
    args = cli.parse_args()
    stats = Counter()
    last = 0
    for line in (
        Path(__file__)
        .with_name("ternary_quadratic_factors.txt")
        .read_text()
        .splitlines()
    ):
        if line.strip() and not line.startswith("#"):
            a, *fields = line.split()
            factors = [
                list(map(int, field.split(","))) if field != "-" else []
                for field in fields
            ]
            if (
                int(a) != last + 1
                or len(fields) != 3
                or any(len(g) != len(set(g)) for g in factors)
            ):
                raise ValueError(("invalid certificate line", line))
            last = int(a)
            verify(last, factors, stats)
    if not last:
        raise ValueError("empty certificate")
    print(f"Independently verified indices 1 through {last}: {dict(stats)}")
    if args.outer_limit is not None:
        verify_outer_limit(args.outer_limit, last)
