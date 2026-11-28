"""Independent quadratic-block checks using linear congruences and integer binomials."""

from argparse import ArgumentParser
from collections import Counter
from functools import lru_cache
from itertools import combinations_with_replacement, product
from math import comb, gcd, isqrt, lcm, prod

from sympy import divisors, factorint


@lru_cache(maxsize=65536)
def factors(n):
    return {int(q): int(q**e) for q, e in factorint(n).items()}


def value(f, x):
    c, b, a = f
    return c + x * (b + a * x)


def all_shapes(b, c):
    # Label the shifted coordinate last, assign the top term, then distribute
    # the two lower coefficients. Exchange the unshifted slots only.
    for top in (1, 2):
        for a0 in range(b + 1):
            for a1 in range(b - a0 + 1):
                for c0 in range(c + 1):
                    for c1 in range(c - c0 + 1):
                        row = (
                            (c0, a0, 0),
                            (c1, a1, int(top == 1)),
                            (c - c0 - c1, b - a0 - a1, int(top == 2)),
                        )
                        if min(map(sum, row)) and (top == 1 or row[0] <= row[1]):
                            yield row


def free(parts, p):
    total = sum(parts)
    w = p
    while w <= total:
        if total // w > sum(x // w for x in parts):
            return False
        w *= p
    return True


def low_line(f, b, c):
    d, k, top = f
    return (b - k, c - d) if top else (k, d)


def line_bound(line, b, c):
    alpha, beta = line
    if not alpha:
        return beta
    r = beta**2 - b * alpha * beta + c * alpha**2
    # In the repeated-root case the line is alpha*(X+h).
    return abs(r) if r else alpha**2


def search_limit(row, B):
    beta, alpha, _ = row[0]
    if not alpha:
        constant = beta
    else:
        d, k, top = row[1]
        constant = abs(d * alpha**2 - k * alpha * beta + top * beta**2)
    # X is above beta. One unshifted part is at most (alpha+1)*X.
    bound = (alpha + 1) ** 3 * B**2 // 27
    if constant:
        bound = min(bound, isqrt(isqrt(constant**3 * B**2 // 27)))
    return bound


@lru_cache(maxsize=262144)
def line_roots(line, w):
    alpha, beta = line
    g = gcd(alpha, w)
    if beta % g:
        return ()
    m = w // g
    r = 0 if m == 1 else -(beta // g) * pow(alpha // g, -1, m) % m
    return tuple(range(r, w, m))


def block_residues(lines, B, b, c):
    residues = [0]
    modulus = 1
    for w in factors(B).values():
        roots = set(line_roots(lines[0], w)) | set(line_roots(lines[1], w))
        roots = [r for r in roots if (r * r + b * r + c) % w == 0]
        inverse = pow(modulus, -1, w)
        residues = [
            s + modulus * ((r - s) * inverse % w) for s in residues for r in roots
        ]
        modulus *= w
    return residues


def check_companion(p, M, parts, B, u, row):
    # For each base separately test the entire power orbit with factorial
    # floors through a fixed depth. Exact rational directions give additional
    # fixed divisors at a pair cut, independently of the power residue.
    for q in sorted({2, 3, 5, 7} | set(factors(M * u))):
        if q == p or u * q > 100000:
            continue
        levels = [q]
        while levels[-1] * u * q <= 100000:
            levels.append(levels[-1] * q)
        modulo = u * levels[-1]
        powers = {1}
        z = p % modulo
        while z not in powers:
            powers.add(z)
            z = z * p % modulo
        for z in powers:
            n = 1 + z * M
            if n % u:
                continue
            v = n // u
            if v % q == 0 and u % q == 0:
                continue
            if B % q == 0 or u % q == 0:
                a, b, c = parts
                if (a + b) // gcd(a, b) % q == 0:
                    continue
                actual = (z * a, z * b, z * c + 1)
            elif M % q == 0 or n % q == 0:
                groups = [(row[0], row[1]), (row[1], row[2]), (row[2], row[0]), row]
                if any(sum(group) // gcd(*group) % q == 0 for group in groups):
                    continue
                actual = tuple(v * x for x in row)
            else:
                break
            if all(sum(x // w for x in actual) == n // w for w in levels):
                break
        else:
            return
    raise ArithmeticError(("surviving exponent classes", p, M, parts, B, u, row))


def inspect(p, M, parts, B, counts):
    d = gcd(*parts[:2])
    remainder = M // B
    if d < 3 or 27 * remainder**2 > d**3 or not free(parts, p):
        return
    denominator = M // gcd(*parts)
    if denominator < 364 or any(
        3 * w >= denominator for w in factors(denominator).values()
    ):
        return
    for u in divisors(d):
        if u < 3 or gcd(u, p * M) > 1 or 27 * remainder**2 > u**3:
            continue
        weights = tuple(factors(remainder).values())
        if any(2 * w + 1 > u for w in weights):
            continue
        allowed = {
            x for x in range(1, u - 1) if all(x % w in (0, u % w) for w in weights)
        }
        for x, y in combinations_with_replacement(sorted(allowed), 2):
            z = u - x - y
            if z < y or z not in allowed:
                continue
            if any(sum(t % w != 0 for t in (x, y, z)) != 1 for w in weights):
                continue
            counts["direct companion tests"] += 1
            check_companion(p, M, parts, B, u, (x, y, z))


def check(b, c, counts):
    start = max(2, b + 1, c + 1)
    for X in range(2, start):
        pf = factors(X)
        if len(pf) != 1:
            continue
        (p,) = pf
        M = X**2 + b * X + c
        for x in range(1, M // 2 + 1):
            for y in range(x, M - x):
                z = M - x - y
                if not free((x, y, z), p):
                    continue
                eligible = [w for w in factors(M).values() if not (x % w and y % w)]
                for choices in product(*[(1, w) for w in eligible]):
                    inspect(p, M, (x, y, z), prod(choices), counts)
    for row in all_shapes(b, c):
        counts["polynomial shapes"] += 1
        lines = tuple(low_line(f, b, c) for f in row[:2])
        cap = lcm(*(line_bound(line, b, c) for line in lines))
        for B in divisors(cap):
            limit = search_limit(row, B)
            if limit < start:
                continue
            for r in block_residues(lines, B, b, c):
                for X in range(start + (r - start) % B, limit + 1, B):
                    M = X**2 + b * X + c
                    if gcd(B, M // B) != 1:
                        continue
                    values = tuple(value(f, X) for f in row)
                    if 27 * (M // B) ** 2 > gcd(*values[:2]) ** 3:
                        continue
                    pf = factors(X)
                    if len(pf) == 1:
                        (p,) = pf
                        inspect(p, M, values, B, counts)


def direct_controls():
    counts = Counter()
    weights = set()
    for h, k in ((0, 2), (1, 3), (2, 3), (2, 4)):
        b, c = h + k, h * k
        for p, a in ((2, 3), (3, 2), (5, 1), (7, 1), (11, 1)):
            X = p**a
            if X <= max(b, c):
                continue
            M = (X + h) * (X + k)
            for e in (1, 2, 3):
                Q = p**e
                n = 1 + Q * M
                if n <= 450:
                    weights.add(n)
                for row in all_shapes(b, c):
                    x, y, z = (value(f, X) for f in row)
                    if not free((x, y, z), p):
                        continue
                    actual = (Q * x, Q * y, 1 + Q * z)
                    if n <= 1500:
                        coeff = comb(n, actual[0]) * comb(n - actual[0], actual[1])
                        if coeff % p == 0:
                            raise ArithmeticError(("digit split", p, n, row))
                        counts["exact carry-free coefficients"] += 1
                    u = gcd(n, gcd(x, y))
                    for root in (h, k):
                        if root < 2 or all(value(f, -root) == 0 for f in row):
                            continue
                        groups = [1, 1, 1]
                        for q, w in factors(X + root).items():
                            full = factors(Q * M)[q]
                            residues = [v % full for v in actual]
                            if sorted(residues) == [0, 0, 1]:
                                groups[residues.index(1)] *= w
                        for j in (1, 2, 3):
                            ell = (j * e + a // 2) // a
                            gap = j * e - ell * a
                            plus, minus = p ** max(0, gap), p ** max(0, -gap)
                            beta, alpha, _ = row[0]
                            if alpha:
                                numerator = beta**2 - b * alpha * beta + c * alpha**2
                                D = plus * (-beta) ** ell * numerator**j - (
                                    -1
                                ) ** j * minus * alpha ** (ell + 2 * j)
                                if D % u:
                                    raise ArithmeticError(
                                        ("weight determinant", n, u, D)
                                    )
                            else:
                                if beta % u:
                                    raise ArithmeticError(("constant part", n, u))
                            for i in (0, 1):
                                D = (
                                    plus * (-root) ** ell * value(row[i], -root) ** j
                                    - minus
                                )
                                if D % groups[i]:
                                    raise ArithmeticError(
                                        ("switched determinant", n, groups[i], D)
                                    )
                                if groups[i] > 1:
                                    counts["nontrivial switched divisibilities"] += 1
                            if any(value(f, -root) % groups[2] for f in row):
                                raise ArithmeticError(("same slot", n, groups))
                            counts["determinant systems"] += 1
    for n in sorted(weights):
        endpoint = n * (n - 1)
        classes = set()
        for x in range(1, n // 3 + 1):
            for y in range(x, (n - x) // 2 + 1):
                classes.add(gcd(endpoint, comb(n, x) * comb(n - x, y)))
        if any(gcd(i, j) == 1 for i, j in combinations_with_replacement(classes, 2)):
            raise ArithmeticError(("exact endpoint pair", n))
    counts["weights with exact endpoint-pair check"] = len(weights)
    if not counts["nontrivial switched divisibilities"] or not weights:
        raise ArithmeticError("missing positive controls")
    return counts


def two_forced_roots():
    # Construct the last possible pair of unshifted polynomials *before*
    # imposing carry freedom. Check their coefficient blocks and the full
    # trinomial valuation separately.
    checked = 0
    for p in (2, 3, 5, 7, 11, 13):
        powers = [1]
        while powers[-1] * p <= 64:
            powers.append(powers[-1] * p)
        for d, r, extra in product(powers, repeat=3):
            k = r + extra
            if d in (r, k):
                continue
            X = p
            while X <= max(d + k, d * k):
                X *= p
            for alpha in range(1, extra + 1):
                beta = extra * d - alpha * r
                gamma = extra - alpha
                if beta < 0 or beta + gamma == 0:
                    continue
                if free((alpha, gamma), p) and free((alpha * r, beta), p):
                    raise ArithmeticError(("residual blocks", p, d, r, extra, alpha))
                low = alpha * (X + r)
                high = (X + r) * (X + d)
                shifted = gamma * X + beta
                if free((p * low, p * high, 1 + p * shifted), p):
                    raise ArithmeticError(
                        ("simultaneous root assignment", p, d, r, extra, alpha)
                    )
                checked += 1
    return checked


if __name__ == "__main__":
    cli = ArgumentParser(description=__doc__)
    cli.add_argument("--height", type=int, default=10)
    args = cli.parse_args()
    counts = Counter()
    for b, c in product(range(args.height + 1), repeat=2):
        delta = b * b - 4 * c
        if delta > 0 and isqrt(delta) ** 2 == delta:
            continue
        check(b, c, counts)
        counts["quadratics"] += 1
    print(dict(counts))
    print(dict(direct_controls()))
    print(
        f"Simultaneous root assignments excluded by factorial floors: {two_forced_roots()}"
    )
