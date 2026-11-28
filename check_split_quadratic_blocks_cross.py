"""Independent enumeration using coefficient triples, linear congruences and power orbits."""

from collections import Counter
from functools import lru_cache
from itertools import product
from math import gcd, isqrt, lcm, prod

from sympy import divisors, integer_nthroot, primerange

from check_shifted_monic_blocks_cross import (
    all_shapes,
    block_residues,
    factors,
    free,
    inspect,
    line_bound,
    low_line,
    search_limit,
    value,
)


def visit(p, M, parts, counts, bounds=None):
    if not free(parts, p):
        return
    u = gcd(*parts[:2])
    if bounds:
        u = min(u, bounds[0])
    if 27 * (M // gcd(M, lcm(*parts[:2]))) ** 2 > u**3:
        return
    choices = [
        (1, w) for w in factors(M).values() if any(x % w == 0 for x in parts[:2])
    ]
    for selected in product(*choices):
        B = prod(selected)
        if (not bounds or B <= bounds[1]) and 27 * (M // B) ** 2 <= u**3:
            inspect(p, M, parts, B, counts)


def bounded_B(row, h, k, start, caps, counts):
    b, c = h + k, h * k
    lines = tuple(low_line(f, b, c) for f in row[:2])
    bounds = [line_bound(line, b, c) for line in lines]
    cap = gcd(prod(caps.values()), lcm(*bounds))
    for B in divisors(cap):
        limit = search_limit(row, B)
        if limit < start:
            continue
        for r in block_residues(lines, B, b, c):
            for X in range(start + (r - start) % B, limit + 1, B):
                M = (X + h) * (X + k)
                parts = tuple(value(f, X) for f in row)
                if gcd(M // B, B) > 1 or 27 * (M // B) ** 2 > gcd(*parts[:2]) ** 3:
                    continue
                fs = factors(X)
                if len(fs) == 1:
                    (p,) = fs
                    inspect(p, M, parts, B, counts)


def scale_bound(row, h, k):
    beta, alpha, _ = row[0]
    if not alpha:
        return beta
    if beta in (0, alpha * h, alpha * k):
        return alpha
    c, b, m = row[1]
    det = abs(alpha**2 * c - alpha * beta * b + beta**2 * m)
    if not m and (not b or c in (0, b * h, b * k)):
        return gcd(det, b or c)
    return det


@lru_cache(maxsize=8192)
def primes_upto(upper):
    return tuple(primerange(2, upper + 1))


def orbit_search(row, h, k, start, counts, left):
    beta, alpha, _ = row[0]
    r = beta // alpha
    d = row[1][1] - r
    t = h + k - d
    if r == d + alpha:
        return
    at_u = r == 1
    target = -(d - r) * (t - r) if at_u else alpha * (r - d)
    L = (t - r) * lcm(abs(d - t), alpha)
    C = abs(target)
    if not (
        (d == 1 or at_u)
        and beta == alpha * r
        and row[1] == (r * d, r + d, 1)
        and 1 < C <= 4
        and 1 <= L <= 256
    ):
        raise ArithmeticError(("orbit height", row, h, k))
    bad = factors(C)
    for p in bad:
        if len(bad) == 1 and all(free(tuple(f[i] for f in row), p) for i in range(2)):
            left.add((p, alpha, r, d, t, target))
    for a in range(1, 82):
        # Outside a full resonance g <= 2*C^2*X^(1/3). With the checked
        # heights above either size test gives X < 2^81.
        cap = C**a + 1
        maxX = isqrt(L**2 * cap**3 // 27) if at_u else 2 * (L * C * cap) ** 2 // 27
        for p in primes_upto(integer_nthroot(maxX, a)[0]):
            X = p**a
            if X < start or X >= 2**81:
                continue
            parts = tuple(value(f, X) for f in row)
            if not free(parts, p):
                continue
            # Direct multiplication until the residue returns to one.
            z = 1
            best = 1
            while True:
                g = gcd(target * z - 1, X + 1)
                if g != X + 1 or not (len(bad) == 1 and p in bad):
                    best = max(best, g)
                z = p * z % (X + 1)
                if z == 1:
                    break
            M = (X + h) * (X + k)
            limits = (best, L * (X + d)) if at_u else (X + r, L * C * best)
            if 27 * M**2 <= limits[0] ** 3 * limits[1] ** 2:
                visit(p, M, parts, counts, limits)


def check(h, k, counts, left):
    start = max(h + k, h * k) + 1
    for X in range(2, start):
        pf = factors(X)
        if len(pf) != 1:
            continue
        (p,) = pf
        M = (X + h) * (X + k)
        for x in range(1, M // 2 + 1):
            for y in range(x, M - x):
                visit(p, M, (x, y, M - x - y), counts)
    for row in all_shapes(h + k, h * k):
        counts["shapes"] += 1
        caps = {d: lcm(*(abs(value(f, -d)) for f in row[:2])) for d in (h, k)}
        caps = {d: m for d, m in caps.items() if m}
        if not caps:
            if h > 2:
                raise ValueError(("no uncovered root", row))
            continue
        if len(caps) == 2:
            bounded_B(row, h, k, start, caps, counts)
            continue
        u = scale_bound(row, h, k)
        if not u:
            orbit_search(row, h, k, start, counts, left)
            continue
        (d,) = caps
        # At X+d the full assigned part divides caps[d].
        capT = isqrt(u**3 // 27)
        for B in divisors(caps[d]):
            for T in range(max(1, (start + d + B - 1) // B), capT + 1):
                if gcd(B, T) != 1:
                    continue
                X = T * B - d
                pf = factors(X)
                if len(pf) == 1:
                    (p,) = pf
                    visit(p, (X + h) * (X + k), tuple(value(f, X) for f in row), counts)


def check_gcd_estimates():
    counts = Counter()
    for p in (2, 3, 5, 7, 11):
        for a in range(1, 16):
            X = p**a
            for c in range(-16, 17):
                if c == 0:
                    continue
                rest, nu = abs(c), 0
                while rest % p == 0:
                    rest //= p
                    nu += 1
                for e in range(1, 2 * a + 1):
                    g = gcd(X + 1, c * pow(p, e, X + 1) - 1)
                    if rest == 1:
                        n, v = divmod(e + nu, a)
                        full = v == 0 and (-1) ** n == (1 if c > 0 else -1)
                        if full != (g == X + 1) or (not full and (g - 1) ** 3 > X):
                            raise ArithmeticError(("signed power", p, a, c, e, g))
                        counts["full power" if full else "partial power"] += 1
                    else:
                        if g**3 > 8 * c**6 * X:
                            raise ArithmeticError(("two multiples", p, a, c, e, g))
                        counts["other constants"] += 1
    return counts


if __name__ == "__main__":
    counts, left = Counter(), set()
    for h in range(1, 11):
        for k in range(h + 1, min(10 - h, 10 // h) + 1):
            check(h, k, counts, left)
    if left != {(3, 3, 1, 2, 4, -3)}:
        raise ArithmeticError(left)
    print(dict(counts), sorted(left))
    print("Verified gcd bounds:", dict(check_gcd_estimates()))
