"""Certify weights P*(p**j+d) for P=p**a, j>=a and bounded positive d."""

from argparse import ArgumentParser
from collections import Counter
from functools import lru_cache
from math import comb, gcd, isqrt, lcm, prod

from flint import fmpz

from check_six_bit_rays import divisors, normalized_rows
from prime_certificates import primes_through
from search_even_endpoint import avoids, prime_powers


@lru_cache(maxsize=131072)
def prime_base(n):
    z = fmpz(n)
    if z.is_prime():
        return n
    if z.is_perfect_power():
        factors = prime_powers(n)
        if len(factors) == 1:
            return next(iter(factors))
    return 0


@lru_cache(maxsize=16384)
def roots(modulus, difference):
    return [r for r in range(modulus) if (r*(r+difference)-1) % modulus == 0]


def scales(first, second, total):
    g = gcd(first, second)
    movable = lcm(first, second)
    common = gcd(movable, total)
    while common > 1:
        movable //= common
        common = gcd(movable, total)
    return divisors(g), divisors(movable)


def progression(same, mixed, difference, remainder):
    for a in roots(same, difference % same):
        # P = a mod same, P = -remainder mod mixed.
        yield a + same*((-remainder-a)*pow(same, -1, mixed) % mixed)


INITIAL_LIMIT = 10**9


@lru_cache(maxsize=32768)
def bases(n):
    return tuple(prime_powers(n))


def subproducts(powers):
    values = [1]
    for w in powers:
        values += [j*w for j in values]
    return values


def check_complement(n, u, rest, counts, rows=normalized_rows):
    factors = prime_powers(rest)
    if any(2*w+1 > u for w in factors.values()):
        return
    counts["CRT searches"] += 1
    v_other = n//u
    for parts in rows(u, tuple(factors.values())):
        counts["complementary rows"] += 1
        other = tuple(v_other*z for z in parts)
        if all(avoids(n, other, p) for p in tuple(factors)+bases(v_other)):
            raise ArithmeticError(("candidate requires independent verification", n, u, rest, other))


def square_progressions(max_gap, t, x, y):
    shared, separate = scales(x, y, t)
    for v in shared:
        for same in shared:
            if gcd(v, same) != 1:
                continue
            for mixed in separate:
                s = same*mixed
                ceiling = (s*s*v**3-1)//27
                if ceiling <= isqrt(INITIAL_LIMIT)-max_gap:
                    continue
                first_gap = t+(-pow(t, -1, mixed) % mixed)
                for d in range(first_gap, max_gap+1, mixed):
                    if gcd(v, d) != 1:
                        continue
                    lo = max(d+1, isqrt(INITIAL_LIMIT)+1-d)
                    if ceiling < lo:
                        continue
                    for residue in progression(same, mixed, d, d-t):
                        # P+d is a multiple of v.
                        residue += s*((-d-residue)*pow(s, -1, v) % v)
                        modulus = v*s
                        yield d, v, s, range(lo+(residue-lo) % modulus, ceiling+1, modulus)


def square(max_gap, rows=normalized_rows):
    counts = Counter()
    for t in range(2, max_gap+1):
        for x in range(1, t//2+1):
            y = t-x
            for d, v, s, powers in square_progressions(max_gap, t, x, y):
                for P in powers:
                    Q = P+d
                    n = P*Q
                    if (n <= INITIAL_LIMIT or gcd(Q//v, v) != 1
                            or 27*(n-1)**2 > (P*v)**3*s*s):
                        continue
                    counts["progression points"] += 1
                    p = prime_base(P)
                    if not p or d % p == 0:
                        continue
                    counts["full prime-power assignments"] += 1
                    row = (P*x, P*y, P*(Q-t))
                    rest, residual = divmod(n-1, s)
                    if residual:
                        raise ArithmeticError((row, s, n))
                    if gcd(rest, s) != 1 or not all(avoids(n, row, r) for r in (p,)+bases(v)+bases(s)):
                        continue
                    counts["first row passes"] += 1
                    check_complement(n, P*v, rest, counts, rows)
    print(f"Square case through offset {max_gap}: {dict(counts)}", flush=True)
    return counts


def higher_weight(P, p, H, d, counts):
    Q = H+d
    n = P*Q
    q_factors = prime_powers(Q)
    shifted = prime_powers(n-1)
    counts["weights"] += 1
    for v in subproducts(q_factors.values()):
        if 2*v > d:
            continue
        u = P*v
        for s in subproducts(w for w in shifted.values() if v*w <= d):
            rest = (n-1)//s
            if (27*rest*rest > u**3 or 4*v*v*s > d*d or
                    any(2*w+1 > u for r, w in shifted.items() if rest % r == 0)):
                continue
            counts["factor assignments"] += 1
            possible = [x for x in range(v, d, v) if P*x*(P*x-1) % s == 0]
            for i, x in enumerate(possible):
                for y in possible[i:]:
                    if x+y > d:
                        break
                    row = (P*x, P*y, P*(Q-x-y))
                    if all(avoids(n, row, r) for r in (p,)+bases(v)+bases(s)):
                        counts["first row passes"] += 1
                        check_complement(n, u, rest, counts)


def higher(max_gap):
    counts = Counter()
    bound = max_gap**4//432
    for p in primes_through(int(fmpz(bound).root(3))):
        P = p
        while P*p*p <= bound:
            H = P*p
            while H*H <= bound*P:
                for d in range(1, max_gap+1):
                    n = P*(H+d)
                    if d % p and n > INITIAL_LIMIT and 432*(n-1)**2 <= P**3*d**4:
                        higher_weight(P, p, H, d, counts)
                H *= p
            P *= p
    print(f"Higher powers through offset {max_gap}: {dict(counts)}", flush=True)


def cross_check():
    # Enumerate small prime-free rows from exact coefficients without a size cutoff.
    total = mixed_count = shared_count = jointly_used = 0
    positive_roots = 0
    extra_scale = 0
    for P in range(2, 30):
        p = prime_base(P)
        if not p:
            continue
        for Q in range(P+1, 2*P):
            if Q % p == 0:
                continue
            n, d = P*Q, Q-P
            for t in range(2, d+1):
                for x in range(1, t//2+1):
                    y = t-x
                    a, b = P*x, P*y
                    coefficient = comb(n, a)*comb(n-a, b)
                    if coefficient % p == 0:
                        continue
                    same = mixed = 1
                    for r, w in prime_powers(n-1).items():
                        if coefficient % r:
                            if gcd(x, y) % w == 0:
                                same *= w
                            else:
                                mixed *= w
                    v = prod(w for r, w in prime_powers(Q).items() if coefficient % r)
                    assert gcd(x, y) % v == 0 and gcd(v, same*mixed) == 1
                    extra_scale += v > 1
                    shared, separate = scales(x, y, t)
                    assert same in shared and mixed in separate
                    assert (1+t*(d-t)) % mixed == 0
                    assert P % (same*mixed) in progression(same, mixed, d, d-t)
                    total += 1
                    shared_count += same > 1
                    mixed_count += mixed > 1
                    jointly_used += min(same, mixed) > 1
    # Check the full progression, also when both factors exceed one.
    for d in range(2, 21):
        for t in range(2, d+1):
            for x in range(1, t//2+1):
                shared, separate = scales(x, t-x, t)
                for a in shared:
                    for b in separate:
                        if (1+t*(d-t)) % b:
                            continue
                        found = set(progression(a, b, d, d-t))
                        direct = {P for P in range(a*b) if (P*(P+d)-1) % (a*b) == 0
                                  and (P+d-t) % b == 0}
                        assert found == direct
                        if min(a, b) > 1:
                            positive_roots += len(found)
    print(f"Checked {total} exact rows ({shared_count} with shared factors, "
          f"{mixed_count} with mixed factors, {jointly_used} with both, {extra_scale} with extra scale) and "
          f"{positive_roots} CRT residues with both factors nontrivial.", flush=True)


if __name__ == "__main__":
    cli = ArgumentParser(description=__doc__)
    cli.add_argument("max_gap", nargs="?", type=int, default=1000)
    args = cli.parse_args()
    if args.max_gap < 1:
        cli.error("use a positive difference bound")
    cross_check()
    if 2*args.max_gap**2 > INITIAL_LIMIT:
        cli.error("the initial theorem must cover the ranges with d >= p**j")
    square(args.max_gap)
    higher(args.max_gap)
