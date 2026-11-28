"""Verify the divisor sieve for short cofactor intervals above a prime-power multiple."""

from argparse import ArgumentParser
from functools import lru_cache
from itertools import product
from math import gcd, isqrt, prod

from flint import fmpz


@lru_cache(maxsize=256)
def factors(n):
    fs = tuple((int(p), int(e)) for p, e in fmpz(n).factor())
    if prod(p**e for p, e in fs) != n or any(not fmpz(p).is_prime() for p, _ in fs):
        raise ArithmeticError(("factorization", n))
    return fs


def divisors(n):
    values = [1]
    for p, e in factors(n):
        old = values.copy()
        for i in range(1, e + 1):
            values.extend(d * p**i for d in old)
    return values


def fits(z, d, p):
    while z:
        if z % p > d % p:
            return False
        z, d = z // p, d // p
    return True


def sieve(p, P, L, b, h, H):
    """Return the complete set passing all positions through the modulus PL."""
    if L % P or not 0 < b * H < L or h % P == 0:
        raise ValueError("need P dividing L, 0 < bH < L, and P not dividing h")
    T = L // P
    h0 = h % P
    if 2 * h0 > P:
        h0 -= P
    bound = b * P * H**2 // L
    hits = set()
    layers = 0
    for j in range(-((bound + h0) // P), (bound - h0) // P + 1):
        v = h0 + j * P
        for d in divisors(abs(v)):
            if d > H or gcd(p, d) != 1:
                continue
            t = v // d
            for x in range(pow(d, -1, P), b * P, P):
                if x < 3:
                    continue
                ell = (x * d - 1) // P
                z = ell + T * t
                if 0 <= z <= b * d:
                    layers += 1
                    if fits(z, b * d, p):
                        hits.add((d, x))
    return layers, sorted(hits)


def floor_layers(n, part, p, modulus):
    power = p
    while power <= modulus:
        if n // power != part // power + (n - part) // power:
            return False
        power *= p
    return True


def direct_sieve(p, P, L, b, h, H):
    hits = set()
    for d in range(1, H + 1):
        w = h * L + d
        if w % p == 0:
            continue
        for x in range(3, b * P):
            if (x * w - 1) % P == 0 and floor_layers(1 + b * P * w, x * w, p, P * L):
                hits.add((d, x))
    return sorted(hits)


def small_checks():
    cases = hits = 0
    for p, a in [(2, 2), (2, 3), (2, 4), (3, 1), (3, 2), (5, 1), (7, 1)]:
        P = p**a
        for L, b in product((P, p * P, p * p * P), (1, 2, 3)):
            H = (L - 1) // b
            if not H:
                continue
            for h in range(1, 2 * P + 1):
                if h % P == 0:
                    continue
                _, actual = sieve(p, P, L, b, h, H)
                expected = direct_sieve(p, P, L, b, h, H)
                if actual != expected:
                    raise ArithmeticError((p, P, L, b, h, actual, expected))
                cases += 1
                hits += len(expected)
    print(
        f"Direct factorial layers agree in {cases} small windows with {hits} surviving pairs."
    )


def two_blocks():
    for p, a in [(2, 1), (2, 3), (2, 7), (3, 2), (3, 4), (5, 3), (7, 2)]:
        P = p**a
        count = 0
        for d in range(1, P):
            if d % p == 0:
                continue
            x = pow(d, -1, P)
            if x < 3:
                continue
            for h in range(P):
                w = d + h * P
                count += floor_layers(1 + P * w, x * w, p, P * P)
        D = p * (p + 1) // 2
        expected = (D - 1) * D ** (a - 1) - 2
        if p != 2:
            expected -= (p + 3) // 2 * ((p + 1) // 2) ** (a - 1)
        if count != expected:
            raise ArithmeticError(("two-block count", p, a, count, expected))
    print(
        "Exact two-block formula agrees with direct factorial floors at seven prime powers."
    )


def square_boundary():
    checked = 0
    for p in (2, 3, 5, 7, 11, 13):
        S = p
        while S < 7:
            S *= p
        while S <= 20000:
            P = S * S
            H = S + isqrt(S)
            _, hits = sieve(p, P, P, 1, 1, H)
            if hits != [(S + 1, P - S + 1)]:
                raise ArithmeticError(("square boundary", S, hits))
            for t in range(1, 20):
                k = t * P - S - 1
                x = P - S + 1
                n = 1 + P * ((1 + k * P) * P + S + 1)
                g = gcd(x, n)
                if g != gcd(x, k - 1) or (t * t + t + 7) % g:
                    raise ArithmeticError(("scale norm", S, t))
            checked += 1
            S *= p
    print(f"Checked the square boundary and scale norm at {checked} prime powers.")


def crt_geometry():
    powers = [(factors(q)[0][0], q) for q in range(2, 401) if len(factors(q)) == 1]
    cases = [0, 0]
    for p, P in powers:
        for qi, (q, Q) in enumerate(powers):
            if q == p:
                continue
            for r, R in powers[qi + 1 :]:
                if r in (p, q):
                    continue
                w = Q * R
                L = P
                while p * L <= w:
                    L *= p
                d = w - L
                if not 0 < d < L:
                    continue
                n = 1 + P * w
                x = pow(w, -1, P)
                kq, kr = pow(P * R, -1, Q), pow(P * Q, -1, R)
                row = (x * w, kq * P * R, kr * P * Q)
                if sum(row) != n or sum(z // (P * L) for z in row) != 1:
                    continue
                g = gcd(gcd(row[0], row[1]), row[2])
                for u in divisors(g):
                    if u < 3:
                        continue
                    got_leading = row[0] >= P * L
                    if d * d < (u + w % u) * w + L // P or (
                        got_leading and d * d < 4 * u * u * w
                    ):
                        raise ArithmeticError(("gap bound", n, row, u))
                    cases[got_leading] += 1
    if min(cases) == 0:
        raise ArithmeticError("missing branch coverage")
    print(
        f"Checked the leading-layer inequality on {sum(cases)} CRT scales, "
        f"including {cases[1]} with the inverse part above the leading power."
    )


# The large windows use a high quotient h congruent to 1 modulo P.
CASES = [
    (
        2,
        64,
        1,
        0,
        2**39,
        9008,
        [
            (4294967297, 18446744069414584321),
            (150323855353, 18446744070641717833),
            (244813135853, 18446744073031398885),
        ],
    ),
    (
        2,
        128,
        1,
        0,
        2**69,
        2007,
        [
            (18446744073709551617, 340282366920938463444927863358058659841),
            (387381625547900583929, 340282366920938463455468859971606975049),
        ],
    ),
    (
        3,
        40,
        2,
        1,
        16 * isqrt(3**41),
        803,
        [
            (6973568800, 12157665457313536600),
            (20920706404, 12157665453826752199),
            (60047662342, 24315330913457118361),
            (62762119216, 12157665443366398996),
            (71479080200, 24315330917943770558),
            (76709256820, 24315330898936543396),
            (78452649020, 24315330905561433758),
            (90656394424, 12157665436392830194),
        ],
    ),
]


def verify_large(case):
    p, a, b, shift, H, expected_count, expected = case
    P = p**a
    L = P * p**shift
    layers, hits = sieve(p, P, L, b, 1, H)
    if (layers, hits) != (expected_count, expected):
        raise ArithmeticError((case, layers, hits))
    for d, x in hits:
        w = L + d
        n = 1 + b * P * w
        if not floor_layers(n, x * w, p, P * L):
            raise ArithmeticError(("false low-digit hit", d, x))
        if len(factors(w)) <= 2 and gcd(x, n) >= 3:
            # At h=1 the complement of a part taking the leading PL contains
            # two positive multiples of P times the two cofactor powers.
            if b != 1 or x * w < P * L or ((n - x * w) // P) ** 2 >= 36 * w:
                raise ArithmeticError(("first window not resolved", n, d, x))
    if p == 2:
        S = isqrt(P)
        d, x = hits[0]
        if (d, x) != (S + 1, P - S + 1) or gcd(x, 1 + P * (L + d)) != 1:
            raise ArithmeticError("square-root boundary")
        if factors(x) != (
            ((18446744069414584321, 1),)
            if a == 64
            else ((769, 1), (442499826945303593556473164314770689, 1))
        ):
            raise ArithmeticError("factorization of boundary multiplier")
    print(
        f"P={p}^{a}, b={b}, L/P={L//P}, H={H}: "
        f"{layers} signed layer cases, {len(hits)} low-digit hits, first window certified.",
        flush=True,
    )


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument(
        "--case",
        type=int,
        choices=range(3),
        help="run just one large window (numbered 0 to 2)",
    )
    args = parser.parse_args()
    if args.case is None:
        small_checks()
        two_blocks()
        square_boundary()
        crt_geometry()
        for case in CASES:
            verify_large(case)
    else:
        verify_large(CASES[args.case])
