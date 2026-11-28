"""Independent factor-pair and exact factorial-floor checks at the three large windows."""

from argparse import ArgumentParser
from math import gcd, prod

from flint import fmpz

from check_predecessor_cofactor_windows import CASES


def factor_pairs(n, limit):
    fs = tuple((int(q), int(e)) for q, e in fmpz(n).factor())
    if prod(q**e for q, e in fs) != n or not all(fmpz(q).is_prime() for q, e in fs):
        raise ArithmeticError(("invalid factors", n))

    def visit(i, u):
        if i == len(fs):
            yield u, n // u
            return
        q, e = fs[i]
        for _ in range(e + 1):
            if u > limit:
                break
            yield from visit(i + 1, u)
            u *= q

    return visit(0, 1)


def first_carry(n, k, p, bound):
    place = p
    while place <= bound:
        if n // place - k // place - (n - k) // place:
            return place
        place *= p
    return 0


def exponent_counts():
    checked = 0
    for p, r, step, a, b in [
        (2, 3, 2, 5, 1),
        (2, 5, 2, 5, 1),
        (3, 2, 2, 4, 1),
        (3, 2, 2, 4, 2),
        (5, 2, 4, 3, 1),
        (5, 2, 4, 3, 2),
    ]:
        P = p**a
        v = 0
        diff = r**step - 1
        while diff % p == 0:
            v += 1
            diff //= p
        for M in (1, 7, 11):
            for e in (1, 5):
                start = M * r**e
                for t in range(a, 2 * a - v + 1):
                    length = p**t
                    width = v + t - a
                    mod = p ** (a + width)
                    w = start % mod
                    mult = pow(r, step, mod)
                    seen = set()
                    total = 0
                    for _ in range(length):
                        if w in seen:
                            raise ArithmeticError("short multiplicative period")
                        seen.add(w)
                        n = 1 + b * P * w
                        total += any(
                            first_carry(n, x * w, p, mod) == 0
                            for x in range(pow(w, -1, P), b * P, P)
                            if x >= 3
                        )
                        w = w * mult % mod
                    if w != start % mod or total * (2 * p) ** (t - a) > b * length * (
                        p + 1
                    ) ** (t - a):
                        raise ArithmeticError(
                            ("exponent bound", p, r, t, total, length)
                        )
                    checked += 1
    for a in range(3, 9):
        P = 2**a
        period = P * P // 4
        w = 1
        survivors = 0
        for _ in range(period):
            x = pow(w, -1, P)
            if x >= 3 and not first_carry(1 + P * w, x * w, 2, P * P):
                survivors += 1
            w = 3 * w % (P * P)
        if survivors != 6 * 3 ** (a - 3) - 2 or w != 1:
            raise ArithmeticError(("common exponent classes", a, survivors))
    print(f"Checked {checked} complete exponent intervals and their fractional bounds.")


def independent(case):
    p, a, b, shift, H, _, expected = case
    P = p**a
    L = P * p**shift
    hits = set()
    # For the prefix 1 the signed product is 1 modulo P. Use positive
    # factorizations on each side and construct x from the complementary factor.
    for m in range(b * H * H // L + 2):
        for sign in (1, -1):
            product = m * P + sign
            if product <= 0:
                continue
            for d, v in factor_pairs(product, H):
                if v >= P:
                    continue
                for block in range(b):
                    x = (block * P + v) if sign == 1 else ((block + 1) * P - v)
                    if x < 3:
                        continue
                    w = L + d
                    if not first_carry(1 + b * P * w, x * w, p, P * L):
                        hits.add((d, x))
    if hits != set(expected):
        raise ArithmeticError(("factor-pair reconstruction", p, a, hits))
    for d, x in hits:
        f = x * (L + d) // (P * L)
        allowed = {((j - f) * d) % P for j in range(b + 1) if j != f}
        mod = P * P * L
        for k in sorted(allowed | {0, 1, p, P - 1, P, 2 * P}):
            w = (1 + k * P) * L + d
            n = 1 + b * P * w
            # For k=P or 2P a carry may first appear above mod. The
            # terminating case k=0 is treated separately in the theorem.
            if k and gcd(k, p) > 1:
                if not first_carry(n, x * w, p, n):
                    raise ArithmeticError(("nonunit high parameter", p, k, d))
            elif k and ((k % P in allowed) != (first_carry(n, x * w, p, mod) == 0)):
                raise ArithmeticError(("next-block mismatch", p, k, d))
    if p == 2:
        S = 2 ** (a // 2)
        x = P - S + 1
        d = S + 1
        # Use Euclidean division directly on the original weight.
        for z in (0, 1, 2, 3, 17, P - 1, P, P + 1, x + 1, P * x + 1):
            n = 1 + P * ((1 + z * P) * P + d)
            if gcd(n, x) != gcd(z - 1, x):
                raise ArithmeticError(("scale gcd", a, z))
    print(
        f"Independent factor-pair check: {len(hits)} exact surviving pairs at {p}^{a}.",
        flush=True,
    )


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=int, choices=range(3))
    args = parser.parse_args()
    if args.case is None:
        exponent_counts()
    for case in (CASES if args.case is None else [CASES[args.case]]):
        independent(case)
