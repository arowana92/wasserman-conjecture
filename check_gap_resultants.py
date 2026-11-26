"""Check the common-root classification and the odd-part estimates by Sylvester determinants."""

from itertools import combinations_with_replacement
from math import gcd
from random import Random

from flint import fmpz_mat, fmpz_poly


def resultant(f, h):
    n, m = f.degree(), h.degree()
    matrix = fmpz_mat(n + m, n + m)
    for poly, offset, count in ((f, 0, m), (h, m, n)):
        for i in range(count):
            for j, coeff in enumerate(poly.coeffs()):
                matrix[offset+i, i+j] = coeff
    return abs(int(matrix.det()))


def oddpart(value):
    return value // (value & -value)


def check(n, s, j, a, b, c):
    l, r = divmod(n, a)
    w, v = j - c*l, s - b*l
    clear = -min(0, w, w+v)
    x = fmpz_poly([0, 1])
    f = (1 << r) * x**a + 1
    h = (1 << (clear+w)) * x**c + (1 << (clear+w+v)) * x**(c+b) - (1 << clear)
    common = f.gcd(h)
    t = common.degree()
    if t:
        assert 2 * abs(a*s-b*n) == a
        assert b > 0 and (a // gcd(a, b)) % 8 == 4 and 2*t <= a
    rem_res = oddpart(resultant(f // common, h))
    common_value = oddpart(abs(int(common(1 << l))))
    actual = gcd((1 << n) + 1, (1 << j) * ((1 << s)+1) - 1)
    assert common_value * rem_res % actual == 0
    exponents = (0, a*j-c*n, a*(j+s)-(b+c)*n)
    span = max(exponents) - min(exponents)
    assert rem_res**a <= (3**a * (1 << span))**(a-t)
    assert common_value**a <= 1 << ((n+a)*t)
    return t > 0


def separated_blocks():
    checked = 0
    for d in range(1, 31):
        for u in range(d, 121, 2*d):
            for i, j in combinations_with_replacement(range(max(1, u-2*d), u+2*d), 2):
                positions = sorted([0, u, i, i+d, j, j+d])
                if 2 * (positions[-1] - positions[1]) >= 3*d:
                    continue
                _, a, b, start, mid, end = positions
                width = max(b, end-start)
                if 3 * start >= 5 * width:
                    assert (mid-start, end-start) == (a, b)
                    checked += 1
    assert checked
    print(f'Checked {checked} separated patterns, allowing repeated units')


def main():
    cases = []
    # Both signs of the half-unit displacement and all four primitive eighth roots.
    for n in range(2, 120, 4):
        for b in (1, 3, 5, 7):
            for sign in (-1, 1):
                s = (b*n + 2*sign) // 4
                if s > 0 and s % 4 == 2:
                    for c in (2, 6):
                        cases.append((n, s, c*n//4, 4, b, c))
    random = Random(9307)
    for _ in range(4000):
        n, a = random.randint(2, 150), random.randint(1, 10)
        power = n & -n
        s = power * (2*random.randrange(12)+1)
        j = random.randrange(2*n)
        b, c = (a*s+n//2)//n, (a*j+n//2)//n
        cases.append((n, s, j, a, b, c))
    positive = sum(check(*case) for case in cases)
    assert positive > 0
    print(f'Checked {len(cases)} polynomial pairs with {positive} nonconstant common factors')


if __name__ == '__main__':
    main()
    separated_blocks()
