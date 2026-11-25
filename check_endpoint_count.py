"""Check the exact exponent margins and the two divisor identities in endpoint_count.md."""

from fractions import Fraction
from math import gcd, prod

from search_endpoint import atoms, smallest_factors


def exponents():
    if not (25**1000 < 2**999 * 13**520 * 12**480
            and 75**3000 < 2**1997 * 13**520 * 62**2480):
        raise ArithmeticError("binary entropy margins")
    upper = (1 + Fraction(733, 500))**3 * Fraction(5, 4)**6 * 3**4
    if not (Fraction(733, 500)**5 > Fraction(27, 4) and upper**16 < 2**195):
        raise ArithmeticError("weighted binary count")


def check_row(v, row, shift_powers):
    u = sum(row)
    n = v*u
    t = prod(shift_powers)
    s = (n-1)//t
    if (n-1) != s*t or gcd(s, t) != 1 or gcd(u, v) != 1:
        raise ArithmeticError((n, t))
    bins = sorted((t//gcd(y, t), y) for y in row)
    (a, y1), (b, y2), (c, y3) = bins
    if a*b*c != t:
        raise ArithmeticError(bins)
    z1, z2, z3 = (y//(t//d) for d, y in bins)
    if (b*c*z1, a*c*z2, a*b*z3) != (y1, y2, y3):
        raise ArithmeticError((row, bins))
    d = s*a - v*z1
    f, g = d*b - v*a*z2, d*c - v*a*z3
    target = (v*a)**2*z2*z3 - d
    if min(d, f, g) <= 0 or f*g != target:
        raise ArithmeticError(("factor pair", row, t))
    slope, const, numerator = z1*b + z2*a, z3*a*b, s*a*b
    det = numerator*const - slope
    denominator = slope*c + const
    if det <= 0 or det % denominator or (numerator*c + 1) != v*denominator:
        raise ArithmeticError(("linear quotient", row, t))
    if (u*s)**3 < 27*(n-1)**2:
        raise ArithmeticError(("modulus bound", n, u, s))
    return a > 1


def main():
    exponents()
    limit = 50
    sieve = smallest_factors(limit**2)
    total = three_bins = 0
    for u in range(3, limit+1):
        for v in range(3, limit+1):
            if gcd(u, v) != 1:
                continue
            powers = [w for _, w in atoms(u*v-1, sieve)]
            for a in range(1, u//3+1):
                for b in range(a, (u-a)//2+1):
                    c = u-a-b
                    eligible = [w for w in powers if sum(y % w == 0 for y in (a,b,c)) == 2]
                    if eligible:
                        three_bins += check_row(v, (a,b,c), eligible)
                        total += 1

    # All three exceptional bins occur in this larger CRT row.
    row = (52195783, 3210219, 11702862)
    m = 7*sum(row)-1
    factors = []
    p = 2
    while p*p <= m:
        w = 1
        while m % p == 0:
            w *= p
            m //= p
        if w > 1:
            factors.append(w)
        p += 1
    if m > 1:
        factors.append(m)
    if not check_row(7, row, factors):
        raise ArithmeticError("expected three nonempty bins")
    print(f"Checked {total+1} row assignments, {three_bins+1} with three nonempty bins.")
    print("All four rational exponent comparisons passed.")


if __name__ == "__main__":
    main()
