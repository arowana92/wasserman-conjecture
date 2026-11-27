"""Integer checks for simultaneous_exponent_intervals.md."""

from collections import Counter
from fractions import Fraction as F
from itertools import product
from math import factorial, gcd


def check_constants():
    # Sum the exponential series for exp(1) and dominate its tail geometrically.
    exp_upper = sum((F(1, factorial(i)) for i in range(8)), F(0))
    exp_upper += F(1, factorial(8))*F(9, 8)
    sqrt_lower = F(173, 100)
    assert sqrt_lower**2 < 3
    assert exp_upper*(1+F(2, 101))/(3*sqrt_lower) < F(27, 50)
    assert exp_upper**3 < F(50, 27)**5
    assert F(5, 4)**2 < 2
    assert F(6, 5)-1 > F(1, 6)
    # Bound lambda_3, hence the exponent in the odd-prime estimate.
    upper = F(733, 500)
    assert upper**3 > upper**2+1 and upper**125 < 2**69
    assert F(44, 100)+F(69, 125) == F(992, 1000)


def check_slopes():
    # Enumerate triples directly, then encode where each binary bit went.
    for width in range(1, 7):
        bound = 1 << width
        obtained = set()
        for a in range(bound):
            for b in range(bound-a):
                for c in range(bound-a-b):
                    if a & b or a & c or b & c:
                        continue
                    code = tuple(((a >> j) & 1) + 2*((b >> j) & 1)
                                 + 3*((c >> j) & 1) for j in range(width))
                    assert code not in obtained
                    obtained.add(code)
        assert obtained == set(product(range(4), repeat=width))


def valuation(n, r):
    assert n > 0
    h = 0
    while n % r == 0:
        n //= r
        h += 1
    return h


def check_classes():
    tested = 0
    for r in (3, 5, 7, 11, 13, 101, 127):
        for start in (1, 173, r**9+17):
            for stride in (1, 24, r**6):
                step = stride*(r-1)
                h = valuation(2**(r-1)-1, r) + valuation(stride, r)
                for width in (0, 1, 2):
                    block = r**width
                    common = r**h
                    modulus = common*block
                    # Large primes use windows straddling one block boundary.
                    for count in (max(1, block-1), block+1):
                        first = pow(2, start, modulus)
                        ratio = pow(2, step, modulus)
                        x = first
                        next_digits = Counter()
                        for _ in range(count):
                            assert x % common == first % common
                            next_digits[x//common] += 1
                            x = x*ratio % modulus
                        assert len(next_digits) == min(block, count)
                        assert max(next_digits.values()) == (count+block-1)//block
                        tested += 1
    return tested


def check_row_digits():
    # Direct factorial floors at small binary powers, with the odd test prime
    # dividing two slopes and the weight. Check the shifted row digit by digit.
    tested = carry_free = 0
    for a in range(1, 28):
        for b in range(1, 65):
            if a & b:
                continue
            for r in (3, 5, 7, 11, 13):
                if gcd(a, b) % r:
                    continue
                for c in range(1, 35):
                    if c & (a | b) or c % r == 0:
                        continue
                    for exponent in range(1, 15):
                        scale = 1 << exponent
                        parts = (a*scale, b*scale, c*scale+1)
                        total = sum(parts)
                        if total % r:
                            continue
                        level = r
                        remainder_parts = list(parts)
                        valuation_sum = 0
                        carry = False
                        while level <= total*r:
                            valuation_sum += total//level-sum(z//level for z in parts)
                            carry |= sum(z % r for z in remainder_parts) >= r
                            remainder_parts = [z//r for z in remainder_parts]
                            level *= r
                        assert carry == (valuation_sum > 0)
                        tested += 1
                        carry_free += not carry
    assert 0 < carry_free < tested
    return tested, carry_free


if __name__ == "__main__":
    check_constants()
    check_slopes()
    progressions = check_classes()
    rows, survivors = check_row_digits()
    print("Exact branch constants and binary shape counts through six bits passed.")
    print(f"{progressions:,} residue windows and {rows:,} independent row checks "
          f"({survivors:,} carry-free).")
