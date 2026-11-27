"""Check simultaneous one-digit weight residues using integer relations and direct counts."""

from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import gcd, prod

from sympy import Matrix, ZZ
from sympy.polys.matrices.normalforms import smith_normal_decomp


@lru_cache(maxsize=None)
def decomposition(rows):
    matrix = Matrix(rows)
    diagonal, left, right = (x.to_Matrix() for x in smith_normal_decomp(matrix.to_DM().convert_to(ZZ)))
    assert diagonal == left*matrix*right and abs(left.det()) == abs(right.det()) == 1
    divisors = [abs(int(diagonal[i,i])) for i in range(matrix.rank())]
    return divisors, left


def compatible(coefficients, rows):
    divisors, left = decomposition(rows)
    targets = [Fraction(-1,c) for c in coefficients]
    return all(prod(t**int(k) for t,k in zip(targets, left.row(i))) == 1
               for i in range(len(divisors), len(rows)))


@lru_cache(maxsize=None)
def masks(rows):
    matrix = Matrix(rows)
    basis_rows = matrix.T.rref()[1]
    basis = matrix[list(basis_rows), :]
    columns = list(basis.rref()[1])
    interpolate = matrix[:, columns]*basis[:, columns].inv()
    found = []
    for values in product((0,1), repeat=len(columns)):
        if any(values):
            target = interpolate*Matrix(values)
            if all(x in (0,1) for x in target):
                found.append(tuple(int(x) for x in target))
    return found


def ordered_blocks(rows):
    # Removing any available block preserves a valid partition on the remaining rows.
    remaining = tuple(range(len(rows)))
    blocks = []
    while remaining:
        available = masks(tuple(rows[i] for i in remaining))
        if not available:
            return None
        choice = available[0]
        blocks.append(tuple(i for i, yes in zip(remaining, choice) if yes))
        remaining = tuple(i for i, yes in zip(remaining, choice) if not yes)
    return blocks


def prescribed_levels(rows, heights):
    for h in sorted(set(heights)):
        active = [i for i, level in enumerate(heights) if level >= h]
        wanted = tuple(int(heights[i] == h) for i in active)
        if wanted not in masks(tuple(rows[i] for i in active)):
            return False
    return True


def density(coefficients, rows, heights, p):
    """Apply the fixed-height formula away from the finite exceptional primes."""
    if not compatible(coefficients, rows) or not prescribed_levels(rows, heights):
        return Fraction(0)
    divisors, left = decomposition(rows)
    targets = [(-pow(c, -1, p)) % p for c in coefficients]
    kernel = 1
    for i in range(len(rows)):
        target = prod(pow(t, int(j), p) for t,j in zip(targets, left.row(i))) % p
        if i >= len(divisors):
            if target != 1:
                return Fraction(0)
        else:
            g = gcd(divisors[i], p-1)
            if pow(target, (p-1)//g, p) != 1:
                return Fraction(0)
            kernel *= g
    exponent = sum(Matrix([row for row,h in zip(rows, heights) if h >= j]).rank()
                   for j in range(1, max(heights)+1))
    return Fraction(kernel, (p-1)**len(divisors)*p**exponent)


def count_directly(coefficients, rows, heights, p):
    modulus = p**(1+max(heights))
    units = [x for x in range(1, modulus) if x % p]
    success = 0
    for values in product(units, repeat=len(rows[0])):
        for c, row, h in zip(coefficients, rows, heights):
            w = prod(x**a for x,a in zip(values, row))
            if (1+c*w) % (p**(h+1)) != p**h:
                break
        else:
            success += 1
    return Fraction(success, len(units)**len(rows[0]))


def check_densities():
    systems = [
        ((1,1,1), ((1,0),(0,1),(1,2))),
        ((4,6,9), ((2,2),(3,4),(4,6))),
        ((1,1), ((1,),(3,))),
        ((1,1,1), ((1,0),(0,1),(1,1))),
        ((2,3,5), ((2,0),(1,1),(0,2))),
    ]
    cases = 0
    for coefficients, rows in systems:
        for p in (5,7):
            if any(c % p == 0 for c in coefficients):
                continue
            for heights in product((1,2), repeat=len(rows)):
                expected = density(coefficients, rows, heights, p)
                actual = count_directly(coefficients, rows, heights, p)
                assert expected == actual, (coefficients, rows, heights, p, expected, actual)
                cases += 1
    print(f"Checked {cases} specified level patterns against all unit tuples.")


def chains():
    for m in range(2,8):
        rows = [tuple(int(j == k) for j in range(m)) for k in range(m)]
        for k in range(m-1):
            rows.append(tuple(int(j == k)+2*int(j == k+1) for j in range(m)))
        rows = tuple(rows)
        levels = tuple(range(1, m+1))+tuple(range(1, m))
        assert compatible((1,)*len(rows), rows) and prescribed_levels(rows, levels)
        blocks = ordered_blocks(rows)
        assert blocks is not None and len(blocks) == m
        assert blocks == [tuple(j for j,h in enumerate(levels) if h == k) for k in range(1,m+1)]
    print("Checked the forced ordering through seven distinct valuation levels.")


def equal_heights():
    rows = ((2,2),(3,4),(4,6))
    assert masks(rows) == [(1,1,1)]
    assert density((4,6,9), rows, (1,1,1), 13) == Fraction(2, 12**2*13**2)
    x, y = 89*193, 89**2*193**3
    assert (1+4*x*x, 1+6*x*y, 1+9*y*y) == (
        1_180_197_317, 5_868_816_171_461_815, 29_184_105_731_699_278_143_082,
    )
    allowed = []
    for a in range(12):
        for b in range(156):
            s, q = pow(89, a, 169), pow(193, b, 169)
            x, y = s*q, s*s*q**3
            if all(n % 169 == 13 for n in (1+4*x*x, 1+6*x*y, 1+9*y*y)):
                allowed.append((a,b))
    assert pow(89,12,169) == pow(193,156,169) == 1
    assert allowed == [(1,1),(1,79)]
    print("Prime 13 certifies all three weights at a=1 mod 12 and b=1 mod 78.")


def concrete():
    from flint import fmpz
    from check_scaled_row_blocks import certify_prime

    p = 251
    bases = []
    for c, h in [(5,1),(7,2)]:
        m = p**(h+1)
        base = (p**h-1)*pow(c, -1, m) % m
        while not fmpz(base).is_prime():
            base += m
        certify_prime(base)
        bases.append(base)
    s, q = bases
    rows = ((1,0),(0,1),(1,2))
    assert prescribed_levels(rows, (1,2,1)) and not prescribed_levels(rows, (1,1,1))
    assert density((5,7,245), rows, (1,2,1), p) == Fraction(1, 250**2*251**3)
    period = p*p*(p-1)
    for u,v in [(1,1),(1+period,1),(1,1+period),(1+period,1+2*period)]:
        x, y = pow(s,u,p**3), pow(q,v,p**3)
        assert (1+5*x) % p**2 == p and (1+7*y) % p**3 == p*p
        assert (1+245*x*y*y) % p**2 == p
    print(f"At p=251 the pair s={s}, q={q} works at all exponents congruent to one modulo {period}.")


if __name__ == "__main__":
    check_densities()
    chains()
    equal_heights()
    concrete()
