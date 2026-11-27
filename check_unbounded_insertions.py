"""Finite checks for unequal widths and unbounded binary gap lengths."""

from collections import Counter
from fractions import Fraction
from itertools import product
from math import gcd

from check_binary_gap_bounds import entropy_below


def parameters():
    counts = Counter()
    for u in range(2, 8):
        for v in range(1, u):
            if gcd(u, v) > 1:
                continue
            for num, den in ((1, 1), (1, 3), (3, 1), (3, 5), (5, 7)):
                t = Fraction(num, den)
                for j, k, sign, delta in product(range(-2, 3), range(-2, 3), (-1, 1), (0, 1)):
                    d = den**u * 4
                    m = den**v * 4
                    x = -Fraction(2)**j * t**u
                    y = sign * Fraction(2)**k * t**v
                    c = delta - d*x - m*y
                    if c.denominator != 1 or c < 2:
                        continue
                    c = int(c)
                    from_den = ((c-delta)*den**u + m*sign*Fraction(2)**k*num**v*den**(u-v)) / (Fraction(2)**j*num**u)
                    from_num = delta + d*Fraction(2)**j*num**u/den**u - m*sign*Fraction(2)**k*num**v/den**v
                    assert from_den == d and from_num == c
                    assert (c-delta) % num == d % den == 0
                    counts["rational points"] += 1
                    counts["negative middle sign" if sign < 0 else "positive middle sign"] += 1
    return counts


def convex_arcs():
    count = 0
    for side in (50, 300, 1500):
        radius = next(r for r in range(1, side + 1) if r**3 >= side)
        for degree in (2, 3, 4, 6):
            for denom in (1, 4, 27, 50, 125, 343):
                # Translated, steep and very flat rational convex graphs.
                offset = side // 3
                points = []
                for x in range(side + 1):
                    y, rem = divmod((x+offset)**degree - offset**degree, denom)
                    if y > side:
                        break
                    if not rem:
                        points.append((x, y))
                steps = [(x2-x1, y2-y1) for (x1,y1),(x2,y2) in zip(points, points[1:])]
                slopes = [Fraction(dy, dx) for dx,dy in steps]
                assert all(a < b for a,b in zip(slopes, slopes[1:]))
                small = sum(dx+dy <= radius for dx,dy in steps)
                assert small <= radius**2
                assert (len(steps) - small)*radius <= 2*side
                count += 1
    return count


if __name__ == "__main__":
    assert entropy_below(3,5,39,40) and entropy_below(1,5,29,40)
    assert Fraction(29,40) + Fraction(9,40) < Fraction(39,40)
    assert Fraction(8,5) - Fraction(2,3)*Fraction(39,40) < Fraction(39,40)
    # The affine-triple exponent is below 5/9, by x^3=x^2+1.
    assert 47**3 > 47**2 * 32 + 32**3 and 47**9 < 2**5 * 32**9
    assert 1 + Fraction(3,32) - (1 - Fraction(5,9))/3 < Fraction(39,40)
    print("Polynomial parametrizations:", dict(parameters()))
    print("Convex graph counts:", convex_arcs())
