"""Exact checks of the two placements and of the resulting common-divisor bound."""
from fractions import Fraction
from itertools import combinations, combinations_with_replacement
from math import gcd


def placement(left, right, shift):
    right = tuple(x+shift for x in right)
    points = sorted(set(left+right))
    gap, high = max((b-a, b) for a, b in zip(points, points[1:]))
    return gap, high, right


def check_pair(left, right, bases):
    height = left[-1]
    if right[-1] > height:
        return check_pair(right, left, bases)
    cases = [placement(left, right, s) for s in (0, height-right[-1])]
    gap, boundary, placed = max(cases, key=lambda item: item[0])
    if 3*right[-1] > 2*height:
        assert 3*gap >= height
    for p in bases:
        x, y = (sum(p**a for a in row) for row in (left, right))
        common = gcd(x, y)
        assert common**3 <= 4**3 * p**(2*height)
        if all(row[0] < boundary <= row[-1] for row in (left, placed)):
            def layers(row):
                return (sum(p**a for a in row if a < boundary),
                        sum(p**(a-boundary) for a in row if a >= boundary))
            a, b = layers(left)
            c, d = layers(placed)
            det = abs(a*d-b*c)
            assert det and det % common == 0 and det <= 4*p**(height-gap)


def constants():
    h = 10**16
    assert h >= 200000**3 and h == 10000**4
    assert Fraction(1, 200)-Fraction(6, 200000)-Fraction(2, h) > Fraction(1, 210)
    # If U and the four-term offset factor are small, a matching
    # factor of degree greater than 11H/35 is compulsory.
    lower = (Fraction(2, 3)-Fraction(3, 200))/2
    assert lower-Fraction(8, 10000)-Fraction(1, 200)-Fraction(2, h) > Fraction(11, 35)
    assert Fraction(49, 11) < 5


if __name__ == '__main__':
    constants()
    shapes = [(0, a, b) for a, b in combinations_with_replacement(range(25), 2) if b]
    exact = binary = 0
    for left, right in combinations(shapes, 2):
        bases = [3, 5, 7]
        if len(set(left)) == len(set(right)) == 3:
            bases.append(2)
            binary += 1
        check_pair(left, right, bases)
        exact += len(bases)
    print(f"{exact} exact gcd bounds, including {binary} binary pairs")

    # Exhaust the geometry separately, including coincident positions.
    geometric = 0
    for h in range(1, 53):
        for a in range(h+1):
            for d in range(2*h//3+1, h+1):
                for c in range(d+1):
                    x, y = (0, a, h), (0, c, d)
                    assert max(placement(x, y, s)[0] for s in (0, h-d))*3 >= h
                    geometric += 1
    print(f"{geometric} normalized endpoint placements")
