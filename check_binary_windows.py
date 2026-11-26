"""Exact constants and finite CRT checks for binary_windows.md (Python 3.9+)."""

from fractions import Fraction as F
from itertools import product


def ln_bounds(q, terms=48):
    """Rational enclosure from the power series for 2 atanh((q-1)/(q+1))."""
    if q < 1:
        lower, upper = ln_bounds(1 / q, terms)
        return -upper, -lower
    t = (q - 1) / (q + 1)
    power = t
    total = F(0)
    for j in range(terms):
        total += 2 * power / (2 * j + 1)
        power *= t * t
    return total, total + 2 * power / ((2 * terms + 1) * (1 - t * t))


LN2_LO, LN2_HI = ln_bounds(F(2))


def log2_upper(q):
    assert q >= 1
    return ln_bounds(q)[1] / LN2_LO


def tangent(t, p):
    """Affine majorant of binary entropy at t, tangent at rational p."""
    return t * log2_upper(1 / p) + (1 - t) * log2_upper(1 / (1 - p))


def constants():
    kappa, lam = F(773, 1000), F(79, 125)
    r0, s0 = F(62256, 100000), F(67094, 100000)
    b = tangent(F(1, 3), F(1, 3))
    upper = 0
    for r, s in product((F(0), F(1)), repeat=2):
        edges = ((2 - lam) * tangent(r, r0)
                 + lam * (b * r + tangent(r / 3, r0 / 3)))
        center = (1 - lam) * tangent(s, s0) + lam * tangent(s / 3, s0 / 3)
        upper = max(upper, (1 - kappa) * edges + (2 * kappa - 1) * center)
    assert upper < F(9329, 10000) < F(14, 15)
    assert -1 + 5 * kappa / 2 < F(9329, 10000)

    theta, ell = F(19, 25), F(7039, 10000)
    d_high, d_low = F(1, 1000), F(1, 5000)
    c = tangent(F(1, 6), F(1, 6))
    cost = (theta - ell + b * (1 - theta + (F(1, 2) + d_high) * (theta - ell))
            + ell * (c + d_low * log2_upper(F(5)) / 3))
    delta = F(1, 10**8)
    assert cost < F(75997, 100000) < theta - delta
    assert -1 + 5 * ell / 2 < theta - delta
    for length, deviation in ((theta - ell, d_high), (ell, d_low)):
        assert 2 * length * deviation**2 / LN2_HI > delta
    print("Global entropy exponent < 0.9329, local typical exponent < 0.75997")
    print("Both local binomial tails save more than 10^-8 in the exponent")


def check_windows():
    cases = 0
    for length in range(3, 13):
        for e in range(1, length):
            digits = length - e
            for h in range(1, digits + 1):
                width = 1 << (digits - h)
                for odd_modulus in range(width | 1, width + 8, 2):
                    high_keys, low_keys = set(), set()
                    for word in range(1 << digits):
                        high = (word % odd_modulus, word >> (digits - h))
                        low = (word % odd_modulus, word % (1 << h))
                        assert high not in high_keys and low not in low_keys
                        high_keys.add(high)
                        low_keys.add(low)
                        cases += 1
    print(f"Both ends identify the word in {cases:,} finite CRT cases")


if __name__ == "__main__":
    constants()
    check_windows()
