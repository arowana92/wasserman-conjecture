"""Check the two-row minor divisibility, equality shapes, and large weights."""

from itertools import product
from math import comb, gcd, isqrt

from check_scaled_row_blocks import certify_prime


def row_data(n, target):
    for a in range(1, n // 3 + 1):
        for b in range(a, (n - a) // 2 + 1):
            row = (a, b, n - a - b)
            coefficient = comb(n, a) * comb(n - a, b)
            minors = [target[j] * row[k] - target[k] * row[j]
                      for j, k in ((1, 2), (2, 0), (0, 1))]
            yield row, coefficient, minors


def exact_coefficient_check(limit=32):
    count = nonzero = reduced = 0
    for n in range(3, limit + 1):
        m = n - 1
        modulus = n * n * m
        for r, u in (((4, 1, 2), (2, 0, 3)), ((0, 0, 2), (1, 3, 1)),
                     ((1, 1, 1), (1, 2, 4)), ((0, 1, 1), (2, 1, 0))):
            st = sum(r) * sum(u)
            for (a, ka, left), (b, kb, right) in product(row_data(n, r), row_data(n, u)):
                common = gcd(ka, kb)
                d, j = gcd(common, n), gcd(common, m)
                multiplier = d * d * j
                for i, l in product(range(3), repeat=2):
                    cell = st * a[i] * b[l] - n * m * r[i] * u[l]
                    value = left[i] * right[l] * cell
                    assert multiplier * value % modulus == 0
                    count += 1
                    nonzero += value != 0
                    reduced += value % modulus != 0
                    # Group the first row at its cut a[i].
                    other_cell = sum(u) * (n - a[i]) * b[l] - n * m * u[l]
                    assert multiplier * a[i] * right[l] * other_cell % modulus == 0
    assert 0 < reduced < nonzero < count
    print(f'{count:,} cells, {nonzero:,} nonzero products, '
          f'{reduced:,} requiring a common endpoint factor. Grouped cuts also passed.')


def window(n, s, t, x, y):
    return x * y * (n * (s * t + s * y + t * x) + x * y) < n * n * (n - 1)


def equality_examples():
    # The two rows pass every lowest assigned layer; higher binary carries
    # exclude them as counterexamples. Test the nonzero equality case itself.
    r, u = (23, 5, 7), (119, 11, 13)
    s, t = sum(r), sum(u)
    passed = 0
    for j in (0, 1, 10**10, 10**20, 10**30, 10**40):
        v = 11972031 + 75150075 * j
        scale = 2 * v - 1
        n, m = v * scale, v * scale - 1
        h, k = (v - 1) // s, (scale + 2) // t
        a = tuple(scale * (h * z + (i == 0)) for i, z in enumerate(r))
        b = tuple(v * (k * z - 2 * (i == 0)) for i, z in enumerate(u))
        assert sum(a) == sum(b) == n
        x = max(abs(s * entry - z * n) for entry, z in zip(a, r))
        y = max(abs(t * entry - z * n) for entry, z in zip(b, u))
        if window(n, s, t, x, y):
            passed += 1
            for i, l in product((1, 2), repeat=2):
                assert s * t * a[i] * b[l] == n * m * r[i] * u[l]
            d, f = scale, -2 * v
            divisor = d + f - 1
            assert divisor and d * (1 - d) % divisor == 0
            assert (n - d) * (n - f) == n * m and d * f == n * divisor
    assert passed >= 3
    print(f'{passed} nonexact progression pairs inside the strict window.')


def large_weights():
    p = 5**168
    q = 42211883872453072075928572096340406768795975134335019729909
    certify_prime(5)
    certify_prime(q)
    length = 10**18
    count = length // 5
    lo, hi = (q * (3 * q + offset) + 1 for offset in (length, 2 * length))
    h, k = hi - 2 * p, hi - 3 * q*q
    assert len(str(lo)) == len(str(hi)) == 118
    assert min(lo - 2 * p, lo - 3 * q*q)**20 > hi**13
    assert h < p and k < q*q and 6 * (h + k) < lo
    x, y = 2 * h, 3 * k
    assert x * y * (hi * (6 + 2 * y + 3 * x) + x * y) < lo**2 * (lo-1)
    assert (2 * h + 2)**3 * 3**6 < (lo-1)**2
    assert length % 5 == 0 and (3*q + length) // 5 > 2**193
    assert count < 2**58 and 2 * length <= 2**61 and 2 * 2**61 < 24 * count
    bound = count - (82 * count - 1) // 193
    assert bound > 110_000_000_000_000_000
    m = 1108391398016982148
    assert length <= m < 2 * length and (3*q + m) % 5**27 == 0 and m < 5**27
    print(f'At least {bound:,} certified weights between the two interval endpoints.')

    q = 7275028100510113835355475391056051625808232326939586246027
    m, e = 4486974951398526273, 101
    certify_prime(q)
    n = q * (e * q + m) + 1
    h, k = n - 2 * p, n - e * q*q
    assert len(str(n)) == 118 and 0 < h < p and 0 < k < q*q
    assert 2 * e * (h + k) < n and window(n, 1, e, h, e * k)
    assert gcd(e, n * (n-1)) == 1 and (n-1) % (q * 5**27) == 0
    assert min(h, k) <= q * 5**27 and min(h, k)**25 > n**16
    print('The separate quotient-101 certificate passed.')


def check_counting_lemma():
    # Check the divisibility before taking logarithms on overlapping intervals.
    from math import factorial, lcm, prod

    def max_power(n):
        largest = 1
        for p in range(2, isqrt(n) + 1):
            power = 1
            while n % p == 0:
                n //= p
                power *= p
            largest = max(largest, power)
        return max(n, largest)

    powers = [0] + [max_power(n) for n in range(1, 601)]
    for length in range(1, 15):
        for cutoff in range(1, 45):
            modulus = factorial(length) * lcm(*range(1, cutoff + 1))
            for start in range(101, 575, 7):
                small = (v for v in range(start, start + length) if powers[v] <= cutoff)
                assert modulus % prod(small) == 0


if __name__ == '__main__':
    exact_coefficient_check()
    equality_examples()
    check_counting_lemma()
    large_weights()
