"""Check tower stabilization and the two binary row transformations independently."""

from itertools import combinations
from math import comb, gcd


def v2(n):
    assert n != 0
    n = abs(n)
    return (n & -n).bit_length() - 1


def tower_mod(x, h, bits):
    """Recursively use the exponent of the unit group of Z/2^bits."""
    if bits == 0:
        return 0
    if h == 1:
        return x % (1 << bits)
    exponent_bits = max(0, bits - 2) if bits != 2 else 1
    exponent = tower_mod(x, h - 1, exponent_bits)
    return pow(x, exponent, 1 << bits)


def fixed_moduli():
    evaluations = distances = levels = 0
    precision = 13
    mask = (1 << precision) - 1
    inputs = range(3, 96, 2)
    for b in (1, 7, 101):
        values = {}
        for x in inputs:
            # At this precision reducing an exponent modulo the original modulus
            # is also valid. This iteration does not recurse to smaller moduli.
            bits = 70
            modulus = 1 << bits
            t = x
            previous = None
            for h in range(1, 19):
                assert t == tower_mod(x, h, bits)
                current = pow(t, b, modulus)
                if previous is not None:
                    shared = min(bits, 2 * (h - 1) + 1)
                    assert (current - previous) % (1 << shared) == 0
                values[h, 0, x] = current & mask
                z = ((current * current - 1) // 8) & mask
                for j in range(1, 38):
                    R = 1 << (j + 2)
                    reference = (pow(current, 1 << j, R << precision) - 1) // R
                    assert reference == z
                    if j <= precision + 3:
                        values[h, j, x] = z
                    if j >= precision:
                        assert z == values[h, precision, x]
                    if h >= precision + 1:
                        assert z == values[precision, min(precision, j), x]
                    z = (z + (1 << (j + 1)) * z * z) & mask
                    levels += 1
                t = pow(x, t, modulus)
                previous = current
                evaluations += 1
        for x, y in combinations(inputs, 2):
            for h in range(1, 19):
                assert v2(values[h, 0, x] - values[h, 0, y]) == v2(x - y)
                distances += 1
                if (x - y) % 4:
                    continue
                for j in (1, 3, 9, precision + 3):
                    assert v2(values[h, j, x] - values[h, j, y]) == v2(x - y) - 2
                    distances += 1
    assert tower_mod(3, 4, 137) == pow(3, 3 ** 27, 1 << 137)
    return evaluations, distances, levels


def normalized_rows():
    comparisons = exact = 0
    for U in range(3, 28, 2):
        for d in range(1, 29):
            if gcd(d, U) != 1:
                continue
            for j in (1, 2, 4, 6, 12, 27):
                R = 1 << (j + 2)
                low = (1 + d) * pow(U, -1, R) % R
                u = (1 + d - U * low) // R
                # Generate Z from the integrality condition on V independently.
                z0 = -(1 + d) * pow(R * d, -1, U) % U
                for Z in (z0, z0 + U):
                    Q = 1 + R * Z
                    assert (1 + d * Q) % U == 0
                    V = (1 + d * Q) // U
                    Y = (d * Z + u) // U
                    assert divmod(V, R) == (Y, low)
                    for f in range(1, U - 1):
                        intercept = f * low // R
                        assert f * V // R == f * Y + intercept
                        if R > d * U + U:
                            assert -U < u <= 0
                            assert intercept == -f * u // U
                        comparisons += 1
                    if U * V <= 240 and V > 0:
                        for f in range(1, U - 1):
                            row = (V * f, V, V * (U - 1 - f))
                            n = sum(row)
                            coefficient = comb(n, row[0]) * comb(n - row[0], row[1])
                            if coefficient % 2:
                                high = [part // R for part in row]
                                assert sum(part % R for part in row) == n % R
                                assert all(not (a & b) for a, b in combinations(high, 2))
                            exact += 1
    return comparisons, exact


if __name__ == "__main__":
    print("Tower evaluations, pairwise distances, normalized power checks:", fixed_moduli())
    print("Normalized part comparisons, exact multinomials:", normalized_rows())
