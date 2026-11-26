"""Exact coefficient, cofactor, polynomial, and prime-power certificate checks."""

from fractions import Fraction
from itertools import product
from math import comb, gcd, prod

from check_scaled_row_blocks import certify_prime


def adjugate(matrix):
    # All inputs are symmetric of order three. Use cyclically oriented minors.
    result = []
    for i in range(3):
        j, k = (i + 1) % 3, (i + 2) % 3
        row = []
        for l in range(3):
            m, n = (l + 1) % 3, (l + 2) % 3
            row.append(matrix[j][m] * matrix[k][n] - matrix[j][n] * matrix[k][m])
        result.append(row)
    return result


def quadratic(matrix, v):
    return sum(v[i] * matrix[i][j] * v[j] for i, j in product(range(3), repeat=2))


def symbolic():
    from sympy import diag, expand, Matrix, symbols
    x, y, a, lam = symbols('x y a lam')
    beta = Matrix([x, y, 1-x-y])
    form = lam * diag(*beta) + (a - lam) * beta * beta.T
    assert expand(form.det() - a * lam**2 * prod(beta)) == 0
    v = Matrix(symbols('r s t'))
    e = sum(v)
    claimed = a * lam * sum(beta[(i+1) % 3] * beta[(i+2) % 3] *
                            (v[i] - e * beta[i])**2 for i in range(3))
    claimed += lam**2 * e**2 * prod(beta)
    assert expand((v.T * form.adjugate() * v)[0] - claimed) == 0


def check_matrix(n, a, b, common):
    m = n - 1
    scale = gcd(n, common)**2 * gcd(m, common)
    denominator = n * n * m
    matrix = []
    for i in range(3):
        row = []
        for j in range(3):
            numerator = a * (n-a) * b[i] * b[j]
            if i == j:
                numerator += a * (a-1) * n * b[i]
            value, remainder = divmod(scale * numerator, denominator)
            assert remainder == 0
            row.append(value)
        matrix.append(row)
    cofactors = adjugate(matrix)
    det = sum(matrix[0][j] * cofactors[0][j] for j in range(3))
    assert det * m**2 * n**3 == scale**3 * a**3 * (a-1)**2 * prod(b)
    assert min(cofactors[i][i] for i in range(3)) > 0 and det > 0
    # The cut is the first component, or its complement by changing the labels.
    t, u = 5, 2
    v = (t-u, -u, -u)
    val = quadratic(matrix, v)
    assert val >= 1 and val * denominator == scale * (
        m * a * (t * b[0] - u * n)**2 + a * (a-1) * t*t * b[0] * (n-b[0]))
    for v in ((-1, 3, -2), (1, 0, 0), (2, 1, 5)):
        e = sum(v)
        val = quadratic(cofactors, v)
        discrepancy = sum(b[(i+1) % 3] * b[(i+2) % 3] * (n*v[i] - e*b[i])**2
                          for i in range(3))
        rhs = a*a * (a-1) * (m * discrepancy + n * (a-1) * e*e * prod(b))
        assert val >= 1 and val * m*m * n**4 == scale**2 * rhs
    return scale


def exhaustive(limit=36):
    count = scaled = 0
    for n in range(3, limit + 1):
        cuts = [(a, comb(n, a)) for a in range(2, n)]
        for i in range(1, n//3 + 1):
            for j in range(i, (n-i)//2 + 1):
                b = (i, j, n-i-j)
                coefficient = comb(n, i) * comb(n-i, j)
                for a, cut in cuts:
                    scale = check_matrix(n, a, b, gcd(coefficient, cut))
                    count += 1
                    # The common endpoint multiplier is necessary in some examples.
                    scaled += (a*(n-a)*b[0]*b[1]) % (n*n*(n-1)) != 0
    assert scaled and scaled < count
    print(f'{count:,} exact binomial and trinomial pairs, {scaled:,} with '
          'a nonintegral unscaled off-diagonal entry at positions 1 and 2.')


def uniform_pair_certificate(exponent, q, expected_count):
    certify_prime(q)
    e, P = 10, 2**exponent
    base = 3 * P
    H = 2**(5 * exponent // 12)
    # n and n+1 together lie in this closed interval.
    low, high = base + H, base + 2*H
    h, k = high-base, high-e*q**3
    assert 2 <= h < P and 0 <= low-e*q**3 <= k < q**3
    assert 3 * e * (h + k) < low and h**7 <= 6**6 * (low-1)**4
    term1 = Fraction(h*h*(h-1)*(e*k)**2, 3*(low-1)*low**2)
    term2 = Fraction(h*h*(h-1)**2*e**2, 27*(low-1)**2)
    assert term1 + term2 < 1
    assert H**5 > 27*(high-1)**2
    assert H*(e*(low-e*q**3))**2 > high**2
    divisor = 210 * q
    step = divisor * 143
    residue = divisor * (-pow(divisor, -1, 143) % 143)
    first = low + (residue-low) % step
    count = 1 + (high-1-first)//step
    assert low <= first < high and first % divisor == (first+1) % 143 == 0
    assert count == expected_count
    print(f'At P=2^{exponent}, {count:,} disjoint consecutive pairs pass.')
    return first


def examples():
    uniform_pair_certificate(420, 933053352955195767836684502058363471888193, 1709179)
    q = 12348858302118435931
    n = uniform_pair_certificate(192, q, 3)
    assert n+1 == 18831305206160042291507368269623000702000273884576598888581
    factors = [q, 2, 2, 3, 5, 7, 47, 420757, 937589, 20194081, 74827799, 129591211]
    assert prod(factors) == n and all(p <= q for p in factors)
    for p in set(factors):
        certify_prime(p)
    assert max(p**factors.count(p) for p in set(factors)) == q
    assert 2*q < n+1-3*2**192


if __name__ == '__main__':
    symbolic()
    exhaustive()
    examples()
