"""Check all bounded directions in the translated-triple finiteness proof."""
from collections import Counter
from functools import cache
from itertools import product
from math import gcd

from flint import fmpz_poly
from sympy import Matrix, Rational

O = (0, 0, 0)
X = (1, 0, 0)
Y = (0, 1, 0)
Z = (0, 0, 1)
ZX = (1, 0, 1)
ZY = (0, 1, 1)
GRID = (O, X, Y, Z, ZX, ZY)
SHAPES = (
    ((X, ZY), (ZX, Y, Z)),
    ((X, Z), (ZX, Y, ZY)),
    ((Z, ZX), (X, Y, ZY)),
)


def sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def systems(shape):
    """Return (binomial character, three unit monomials) for each component."""
    x, y = SHAPES[shape]
    yield sub(*x), y                     # U
    if shape != 2:
        yield sub(*x), (O, X, Y)         # A2
    if shape == 0:
        yield Z, (O, *x)                 # C, for either sign of Z
    elif shape == 1:
        yield sub(Y, Z), (O, *x)
    else:
        yield ZY, y


def offsets(shape, p):
    points = set()
    for char, (u, v, w) in systems(shape):
        mat = Matrix([char, sub(v, u), sub(w, u)])
        assert abs(mat.det()) in (1, 2)
        inv = mat.inv()
        for gamma in ((Rational(1, 2), 1) if p == 2 else (Rational(1, 2),)):
            for leading in range(3):
                heights = [gamma if i == leading else 0 for i in range(3)]
                h = inv * Matrix([0, heights[1]-heights[0], heights[2]-heights[0]])
                # Use a common denominator of four for all offsets.
                scaled = tuple(4 * a for a in h)
                assert all(a.is_Integer for a in scaled)
                points.add(tuple(map(int, scaled)))
    return points


def directions(limit):
    for a, b, c in product(range(-limit, limit + 1), repeat=3):
        height = max(0, a, b) - min(0, a, b) + abs(c)
        if 0 < height <= limit and gcd(a, b, c) == 1:
            yield (a, b, c), height


@cache
def terms_to_poly(p, terms):
    coeffs = [0] * (1 + max(n for n, _ in terms))
    for n, k in terms:
        coeffs[n] += p**k
    return fmpz_poly(coeffs)


def specialize(p, support, v, e):
    powers = [(dot(s, v), dot(s, e)) for s in support]
    lo = min(n for n, _ in powers)
    low = min(k for _, k in powers)
    return terms_to_poly(p, tuple(sorted((n-lo, k-low) for n, k in powers)))


def admissible(shape, v, e):
    # Coordinate equalities on the whole line. The original columns are distinct
    # and the translation is nonzero. Two more identities would force U = 1.
    tests = [X, Y, sub(X, Y), Z]
    if shape == 1:
        tests.append(sub(Y, Z))
    if shape == 2:
        tests.append(ZY)
    return all(dot(t, v) or dot(t, e) for t in tests)


def check_line(p, shape, v, e, height):
    x, y = (specialize(p, s, v, e) for s in SHAPES[shape])
    m = specialize(p, GRID, v, e)
    assert m.degree() == height
    u = x.gcd(y)
    s = m.gcd(x*y)
    return 2*height - 3*u.degree() - 2*s.degree()


def main():
    rays = tuple(directions(36))
    counts = Counter()
    bad = []
    duplicate_lines = 0
    for p, shape in product((2, 3), range(3)):
        points = offsets(shape, p)
        for v, height in rays:
            seen = set()
            if height > (36 if shape == 0 else 21):
                continue
            for h in points:
                # A primitive direction makes the parameter's denominator divide 4.
                for step in range(4):
                    nums = tuple(a + step*b for a, b in zip(h, v))
                    if any(n % 4 for n in nums):
                        continue
                    e = tuple(n // 4 for n in nums)
                    # Give each parallel line a unique lattice representative.
                    i = next(i for i in range(3) if v[i])
                    q = e[i] // v[i]
                    e = tuple(a-q*b for a, b in zip(e, v))
                    if e in seen or e == O or not admissible(shape, v, e):
                        continue
                    seen.add(e)
                    margin = check_line(p, shape, v, e, height)
                    counts[p, shape] += 1
                    if 3*margin < height:
                        pairs = ((a, b) for a in GRID for b in GRID if a != b)
                        if p == 2 and any(dot(sub(a, b), v) == dot(sub(a, b), e) == 0
                                          for a, b in pairs):
                            duplicate_lines += 1
                        else:
                            bad.append((p, shape, v, e, margin))
    print(f"{len(rays)} oriented primitive directions")
    for (p, shape), count in sorted(counts.items()):
        print(f"p={p}, type={shape+1}: {count} lines")
    print(f"{sum(counts.values())} oriented line checks, {duplicate_lines} binary repetitions")
    if bad:
        raise ArithmeticError(("line below the one-third margin", bad))
    print("Every admissible line has degree margin at least one third of its height.")


if __name__ == "__main__":
    main()
