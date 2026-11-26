"""Independently evaluate the product-family projections and binary gap constraints."""

from itertools import combinations

from check_geometric_product import GRID, allocations, sectors
from check_polynomial_powers import resultant


def numeric_coefficients(block, v, w, base):
    terms = [(i*v[0]+j*v[1], i*w[0]+j*w[1]) for i, j in block]
    lo, low = min(a for a, b in terms), min(b for a, b in terms)
    values = [0]*(max(a for a, b in terms)-lo+1)
    for i, j in terms:
        values[i-lo] += base**(j-low)
    return tuple(values)


def projected_values():
    tests = 0
    for x, y, z in allocations():
        pairs = [(x, y), (x, tuple(sorted(y+(z,)))), (y, tuple(sorted(x+(z,))))]
        for v, w, fs, _ in sectors((1, 0), (0, 1), pairs):
            for (a, b), h in zip(pairs, fs):
                stripped = []
                for base in (3, 5):
                    f, g = [numeric_coefficients(block, v, w, base) for block in (a, b)]
                    if len(f) == 1:
                        value = f[0]
                    elif len(g) == 1:
                        value = g[0]
                    else:
                        value = resultant(f, g)
                    actual = abs(sum(c*base**i for i, c in enumerate(h)))
                    assert actual and value % actual == 0
                    quotient = value//actual
                    exponent = 0
                    while quotient % base == 0:
                        quotient //= base
                        exponent += 1
                    assert quotient == 1
                    stripped.append(exponent)
                    tests += 1
                assert stripped[0] == stripped[1]
    return tests


def pairings(bits):
    if not bits:
        yield ()
        return
    for j in range(1, len(bits)):
        for rest in pairings(bits[1:j]+bits[j+1:]):
            yield ((bits[0], bits[j]),)+rest


def binary_checks():
    matches = list(pairings(tuple(range(6))))
    distinct = repetitions = 0
    for i in range(1, 65):
        for k in range(1, 129):
            bits = sorted(row*i+col*k for row, col in GRID)
            m = sum(2**j for j in bits)
            if len(set(bits)) == 6:
                for matching in matches:
                    gaps = [bits[v]-bits[u] for u, v in matching]
                    for n in set(gaps):
                        if gaps.count(n) == 2 and max(gaps)-min(gaps) == 1:
                            assert 27*m*m > (2**n+1)**3
                distinct += 1
            else:
                assert k in (i, 2*i)
                digits = tuple(j for j in range(m.bit_length()) if m & (1 << j))
                if len(digits) >= 5:
                    assert k == 2*i and i >= 2
                    for a, b, c, d in combinations(digits, 4):
                        for left, right in (((a,b),(c,d)), ((a,c),(b,d)), ((a,d),(b,c))):
                            n = left[1]-left[0]
                            if n == right[1]-right[0] and 27*m*m <= (2**n+1)**3:
                                assert n == 3*i and right[0]-left[0] == i
                repetitions += 1
    return distinct, repetitions


if __name__ == '__main__':
    print('Independent Sylvester and constant evaluations', projected_values())
    print('Distinct and coincident binary profiles', binary_checks())
