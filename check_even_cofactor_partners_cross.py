"""Independent factorial-floor and small unit-group checks for the surviving rows."""

from fractions import Fraction
from itertools import combinations
from math import gcd

from check_even_cofactor_partners import CASES
from check_shifted_prefixes import floor_cover


def log_classes(target, modulus):
    """For every possible unit base, find its exponents hitting the target."""
    classes = []
    for base in range(1, modulus):
        if gcd(base, modulus) != 1:
            continue
        first = None
        value = 1
        for exponent in range(1, modulus + 1):
            value = value * base % modulus
            if value == target and first is None:
                first = exponent
            if value == 1:
                break
        assert value == 1
        if first is not None:
            classes.append((first, exponent))
    return classes


def verify(digits):
    expected = (Fraction(1, 720), Fraction(9, 15200), Fraction(1, 45))
    for case, fraction in zip(CASES, expected):
        w, P, u, a, b, initial, r, m = case
        D = w * P
        slopes = (P * a, w * b, (w - a) * P - w * b)
        found, steps = floor_cover(slopes, (2, u, P), D*u, initial, 10**digits)
        assert not found
        odd_u = m // ((m & -m) * P)
        sizes = (odd_u // u * (u-1), P-1)
        t = max(s & -s for s in sizes)
        groups = [log_classes(r % mod, mod) for mod in (odd_u, P)]
        binary_fraction = Fraction(1, 8 if D == 24352 else 4)
        for odd in range(1, t, 2):
            count = 0
            for first in groups[0]:
                for second in groups[1]:
                    congruences = (first, second, (odd, t))
                    if all((i-j) % gcd(s1, s2) == 0
                           for (i, s1), (j, s2) in combinations(congruences, 2)):
                        count += 1
            assert binary_fraction * Fraction(count, sizes[0]*sizes[1]) == fraction
        print(D, "floor intervals", steps, "prime fraction", fraction, flush=True)


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--digits", type=int, default=10000)
    verify(parser.parse_args().digits)
