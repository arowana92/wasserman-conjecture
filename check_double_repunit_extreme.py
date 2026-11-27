"""Independent integer check of a length-26 partition with the largest cutoff."""

from math import gcd, isqrt, lcm

PARTS = ((0, 1, 3, 4, 7), (2, 10, 11), (5, 6, 8, 9, 12))


def norm(poly, positions):
    """Determinant of multiplication modulo a monic polynomial, by Bareiss."""
    n = len(poly) - 1
    vector = [int(i in positions) for i in range(n)]
    assert max(positions) < n
    matrix = []
    for _ in range(n):
        matrix.append(vector)
        top = vector[-1]
        vector = [a - top * b for a, b in zip([0] + vector[:-1], poly)]
    previous = 1
    for i in range(n - 1):
        if not matrix[i][i]:
            j = next((j for j in range(i + 1, n) if matrix[j][i]), None)
            if j is None:
                return 0
            matrix[i], matrix[j] = matrix[j], matrix[i]
        pivot = matrix[i][i]
        for j in range(i + 1, n):
            for k in range(i + 1, n):
                numerator = matrix[j][k] * pivot - matrix[j][i] * matrix[i][k]
                assert numerator % previous == 0
                matrix[j][k] = numerator // previous
        previous = pivot
    return abs(matrix[-1][-1])


def group_values(q):
    return tuple(sum(q**j for j in positions) for positions in PARTS)


def second_largest_divisibility(a, b, c):
    return lcm(gcd(a, b), gcd(b, c), gcd(c, a))


def digit_sum(value, p):
    s = 0
    while value:
        value, d = divmod(value, p)
        s += d
    return s


def check():
    # Degree 12 may occur in the inputs, so first reduce x^12 modulo Phi_13.
    phi13 = [1] * 13

    def weight_norm(positions):
        # Multiplication by X^j for a cyclic shift preserves the absolute norm.
        for j in range(13):
            shifted = tuple((i + j) % 13 for i in positions)
            if max(shifted) < 12:
                return norm(phi13, shifted)

    assert [weight_norm(f) for f in PARTS] == [53, 53, 53]
    phi25 = [int(i % 5 == 0) for i in range(21)]
    assert [norm(phi25, f) for f in PARTS] == [63125, 4651, 63125]
    assert 63125 == 5**4 * 101
    assert all(sum(q**i for i in range(0, 21, 5)) % 25 for q in range(25))

    modulus = 53 * 505
    small = [p for p in range(2, 501) if all(p % i for i in range(2, isqrt(p) + 1))]
    tried = carried = largest = 0
    for r in range(modulus):
        gs = group_values(r)
        u = gcd(53, *gs)
        c = gcd(505, sum(r**i for i in range(0, 21, 5)), second_largest_divisibility(*gs))
        limit = ((2**13 + 1) * u) ** 3 * c * c // (27 * 2**39)
        largest = max(largest, limit)
        for q in range(2 + (r - 2) % modulus, limit + 1, modulus):
            tried += 1
            gs = group_values(q)
            scale = q**13 + 1
            rows = tuple(scale * g for g in gs)
            weight = sum(rows)
            upper = gcd(*rows)
            lower = (weight - 1) // gcd(weight - 1, second_largest_divisibility(*gs))
            if 27 * lower**2 > upper**3:
                continue
            carried += 1
            for p in small:
                if upper % p and (weight - 1) % p:
                    continue
                if sum(digit_sum(x, p) for x in rows) == digit_sum(weight, p):
                    continue
                while upper % p == 0:
                    upper //= p
                remainder = weight - 1
                while remainder % p == 0:
                    remainder //= p
                lower = lcm(lower, (weight - 1) // remainder)
                if 27 * lower**2 > upper**3:
                    break
            assert 27 * lower**2 > upper**3, (q, rows)
    print("Largest cutoff", largest, "integer bases tested", tried, "needing carries", carried)


if __name__ == "__main__":
    check()
