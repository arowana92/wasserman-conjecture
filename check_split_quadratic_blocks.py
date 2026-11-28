"""Unrestricted shifts at split monic quadratics of coefficient height ten."""

from collections import Counter
from functools import lru_cache
from math import gcd, isqrt, lcm, prod

from sympy import integer_nthroot

from check_close_powers import prime_base, subproducts
from check_polynomial_powers import divisor_bound, evaluate, power_limit, row_shapes
from check_shifted_monic_blocks import (
    coprime_part,
    inspect_row,
    root_classes,
    small_inputs,
)
from check_six_bit_rays import divisors
from prime_certificates import primes_through
from search_even_endpoint import avoids, prime_powers


def divide_linear(f, d):
    # Divide by X+d when it divides exactly.
    if len(f) < 2 or evaluate(f, -d):
        return None
    out = [f[-1]]
    for c in reversed(f[1:-1]):
        out.append(c - d * out[-1])
    return tuple(reversed(out))


def strip(f, roots):
    for d in roots:
        while (quotient := divide_linear(f, d)) is not None:
            f = quotient
    return f


def numerical(X, p, row, caps, counts):
    values = tuple(evaluate(f, X) for f in row)
    M = sum(values)
    u = coprime_part(gcd(*values[:2]), p * M)
    if u < 3 or not avoids(M, values, p):
        return
    bound = gcd(M, lcm(*values[:2]))
    if 27 * (M // bound) ** 2 > u**3:
        return
    eligible = [
        w for w in prime_powers(M).values() if any(v % w == 0 for v in values[:2])
    ]
    for B in subproducts(eligible):
        if caps is not None and (
            B > caps[1] or 27 * (M // B) ** 2 > min(u, caps[0]) ** 3
        ):
            continue
        if 27 * (M // B) ** 2 <= u**3:
            counts["numerical assignments"] += 1
            inspect_row(p, M, values, B, counts)


def bounded_scale(h, k, row, u, covered, low, counts):
    if u < 3:
        return
    d = min(covered, key=covered.get)
    limit = isqrt(u**3 // 27)
    candidates = {b * t - d for b in divisors(covered[d]) for t in range(1, limit + 1)}
    for X in candidates:
        if X >= low and (p := prime_base(X)):
            numerical(X, p, row, None, counts)


def bounded_predecessor(h, k, row, upper_u, covered, low, counts):
    G = (h * k, h + k, 1)
    # B is a divisor of cap_h*cap_k, with possible shared prime factors.
    for B in divisors(prod(covered.values())):
        high = power_limit(G, 1, upper_u, (B,))
        if high is None:
            raise ArithmeticError(("nonconstant cutoff", row))
        if high < low:
            continue
        powers = tuple(prime_powers(B).values())
        for r in root_classes(h + k, h * k, B):
            if any(evaluate(row[0], r) % w and evaluate(row[1], r) % w for w in powers):
                continue
            for X in range(low + (r - low) % B, high + 1, B):
                M = evaluate(G, X)
                if gcd(B, M // B) > 1:
                    continue
                vals = tuple(evaluate(f, X) for f in row)
                if 27 * (M // B) ** 2 > gcd(*vals[:2]) ** 3:
                    continue
                if p := prime_base(X):
                    inspect_row(p, M, vals, B, counts)


@lru_cache(maxsize=100000)
def max_orbit_gcd(p, a, c, nonfull):
    X = p**a
    return max(
        (
            g
            for t in range(2 * a)
            if (g := gcd(X + 1, c * pow(p, t, X + 1) - 1)) < X + 1 or not nonfull
        ),
        default=1,
    )


def critical(h, k, row, low, counts, residual):
    x, y, z = row
    if len(x) > len(y):
        x, y = y, x
    if len(x) != 2 or len(y) != 3 or x[0] % x[1]:
        raise ArithmeticError(("expected one common linear factor", row))
    alpha, r = x[1], x[0] // x[1]
    d = y[1] - r
    other = h + k - d
    if not (d in (h, k) and y == (r * d, r + d, 1) and 0 < r < other and r != d):
        raise ArithmeticError(("root recovery", row))
    if r == d + alpha:
        counts["primitive pair sum"] += 1
        return
    at_u = r == 1
    if not at_u and d != 1:
        raise ArithmeticError(("no unit root", row))
    c = -(d - r) * (other - r) if at_u else alpha * (r - d)
    C = abs(c)
    if C <= 1:
        raise ArithmeticError(("unit target", row, c))
    L = (other - r) * lcm(alpha, abs(other - d))
    common = 8 * L**2 * C**6
    cutoff = common // 27 if at_u else common**3 // 27**3
    resonant_primes = prime_powers(C)
    if len(resonant_primes) == 1:
        (p,) = resonant_primes
        if all(
            avoids(b, tuple((f[i] if i < len(f) else 0) for f in row), p)
            for i, b in enumerate((h * k, h + k))
        ):
            residual.add((h, k, p, d, r, alpha, c))
    counts["unit-root orbit shapes"] += 1
    for a in range(1, cutoff.bit_length() + 1):
        # The exponent-independent resultant bounds the base even at a=1.
        base_limit = (
            8 * L**2 * C ** (3 * a) // 27 if at_u else 8 * L**2 * C ** (2 * a + 2) // 27
        )
        upper = min(
            integer_nthroot(base_limit, 2 * a if at_u else a)[0],
            integer_nthroot(cutoff, a)[0],
        )
        for p in primes_through(upper):
            X = p**a
            if X < low:
                continue
            nonfull = len(resonant_primes) == 1 and p in resonant_primes
            g = max_orbit_gcd(p, a, c, nonfull)
            u_bound, b_bound = (g, L * (X + d)) if at_u else (X + r, L * C * g)
            if 27 * ((X + h) * (X + k)) ** 2 <= u_bound**3 * b_bound**2:
                counts["orbit specializations"] += 1
                numerical(X, p, (x, y, z), (u_bound, b_bound), counts)


def check(h, k, counts, residual):
    G = (h * k, h + k, 1)
    small_inputs(G, counts)
    low = max(G) + 1
    for shape in row_shapes(G):
        for shifted in range(3):
            if shifted == 1 and shape[0] == shape[1]:
                continue
            unshifted = tuple(f for i, f in enumerate(shape) if i != shifted)
            row = (*unshifted, shape[shifted])
            covered = {}
            for d in (h, k):
                vals = [abs(evaluate(f, -d)) for f in unshifted]
                if all(vals):
                    covered[d] = lcm(*vals)
            counts["shapes"] += 1
            if not covered:
                if h > 2:
                    raise ValueError("both roots require a separate scale test")
                counts["scale at most two"] += 1
                continue
            caps = [
                strip(f, (0, h, k)) for f in (*unshifted, divisor_bound(*unshifted))
            ]
            upper_u = min(caps, key=lambda f: (len(f), sum(map(abs, f))))
            if len(upper_u) == 1:
                bounded_scale(h, k, row, abs(upper_u[0]), covered, low, counts)
            elif len(covered) == 2:
                bounded_predecessor(h, k, row, upper_u, covered, low, counts)
            else:
                critical(h, k, row, low, counts, residual)
    counts["quadratics"] += 1


if __name__ == "__main__":
    counts = Counter()
    residual = set()
    for h in range(1, 11):
        for k in range(h + 1, 11):
            if h * k <= 10 and h + k <= 10:
                check(h, k, counts, residual)
    expected = {(2, 4, 3, 2, 1, 3, -3)}
    if residual != expected:
        raise ArithmeticError(("unexpected resonant shapes", residual))
    print(dict(counts))
    print("Remaining full power orbit:", sorted(residual))
