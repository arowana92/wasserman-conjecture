"""Exact finite checks of the digit-insertion counting argument."""

from collections import Counter, defaultdict
from fractions import Fraction
from itertools import product
from math import factorial, prod

from check_binary_gap_bounds import entropy_below
from check_three_binary_blocks import matrix_rank


def trailing_bit(value):
    return value & -value


def operations():
    checks = Counter()
    widths = (3, 2, 3)
    fields = list(product(*(range(1 << w) for w in widths)))
    for order in ((0, 1, 2), (0, 2, 1), (2, 1, 0)):
        ws = [widths[i] for i in order]
        for reverse in (False, True):
            for flip in (False, True):
                for gap, pattern in ((0, 0), (2, 1), (5, 20)):
                    shifts = (0, ws[0] + gap, ws[0] + gap + ws[1] + 1)
                    fixed = pattern << ws[0]
                    outputs = set()
                    for vals in fields:
                        m = vals[0] + (vals[1] << 3) + (vals[2] << 5)
                        n = fixed
                        for i, pos in enumerate(order):
                            v = vals[pos]
                            if i:
                                if reverse:
                                    v = int(format(v, f"0{ws[i]}b")[::-1], 2)
                                if flip:
                                    v ^= (1 << ws[i]) - 1
                            n += v << shifts[i]
                        assert n not in outputs
                        outputs.add(n)
                        if order[0] == 0 and vals[0] and m % 2 == 0:
                            assert trailing_bit(n) == trailing_bit(m)
                            checks["retained lowest one"] += 1
                        if n >= 3 and n % 2:
                            p = trailing_bit(n-1)
                            t = (n-1)//p
                            assert t & 1 and p*t+1 == n
                            checks["odd projected weight"] += 1
                    checks["injective descriptions"] += 1
    return checks


def field_lattices():
    checks = Counter()
    for widths in ((1, 1, 2, 2), (2, 1, 1, 2, 1), (1, 1, 1, 1, 1, 1)):
        sides = [1 << w for w in widths]
        dim = len(widths)
        points = list(product(*(range(q) for q in sides)))
        for gaps in ((0,)*dim, tuple(i % 3 for i in range(dim))):
            coeff = []
            bit = 0
            for gap, width in zip(gaps, widths):
                bit += gap
                coeff.append(1 << bit)
                bit += width
            def encode(vector):
                return sum(c*z for c,z in zip(coeff, vector))
            possible = set(map(encode, points))
            assert len(possible) == len(points)
            # This also checks every nonzero difference in these disjoint fields.
            checks["noncancelling field differences"] += len(points)*(len(points)-1)
            for mod in (3, 9, 15, 25, 35, 101):
                groups = defaultdict(list)
                for p in points:
                    groups[encode(p) % mod].append(p)
                scale = max(sides)
                def length(x, y):
                    return max(abs(a-b)*scale//w for a,b,w in zip(x,y,sides))
                near = min((length(x,y) for g in groups.values()
                            for j,x in enumerate(g) for y in g[j+1:]), default=scale)
                for group in groups.values():
                    rank = matrix_rank(([z-b for z,b in zip(y,group[0])] for y in group),dim)
                    if rank < dim:
                        # Euclidean packing argument, replacing sqrt(dim) by dim.
                        assert len(group)*near**rank <= (near + 2*dim*scale)**rank
                        checks["lower-dimensional slices"] += 1
                    else:
                        assert len(group)*mod <= 3**dim * factorial(dim)*prod(sides)
                        checks["full-dimensional slices"] += 1
    return checks


def margins():
    assert entropy_below(3,5,39,40) and entropy_below(1,5,29,40)
    for r in range(1, 4097):
        e = Fraction(1, 12*r+32)
        assert Fraction(39,40) < 1-e
        assert Fraction(29,40) + Fraction(1,100) + 3*r*e < 1-e
        assert Fraction(29,40) + Fraction(1,100) + Fraction(1,5) < 1-e
        assert Fraction(8,5) - Fraction(2,3)*(1-e) < 1-e
        assert Fraction(4,5) + Fraction(1,6) < 1-e
        assert e > Fraction(1, 12*r+33)


if __name__ == "__main__":
    margins()
    print("Digit operations:", dict(operations()))
    print("Lattice checks:", dict(field_lattices()))
