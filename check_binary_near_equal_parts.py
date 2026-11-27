"""Check least gaps against exact positive trinomial coefficients."""

from itertools import combinations
from math import comb, factorial


def check_weights(limit):
    factorials = [factorial(i) for i in range(limit + 1)]
    odd_rows = 0
    bit_minima = 0
    for n in range(3, limit + 1):
        actual = {}
        for a in range(1, n // 3 + 1):
            for b in range(a, (n - a) // 2 + 1):
                c = n - a - b
                coefficient = factorials[n] // (factorials[a] * factorials[b] * factorials[c])
                if coefficient % 2 == 0:
                    continue
                odd_rows += 1
                for v, u in combinations((a, b, c), 2):
                    assert v < u
                    t = 1 << (u.bit_length() - 1)
                    gap = u - v
                    delta = 2 * t - n % (2 * t)
                    assert 1 <= delta <= gap
                    actual[t] = min(actual.get(t, gap), gap)
        bits = [1 << i for i in range(n.bit_length()) if n & (1 << i)]
        expected = {}
        if len(bits) >= 3:
            for t in bits[1:]:
                expected[t] = t - n % t + (bits[0] if t == bits[-1] else 0)
            assert min(actual.values()) == bits[1] - bits[0]
        assert actual == expected, (n, actual, expected)
        bit_minima += len(actual)
    print(f"Weights 3 through {limit}: {odd_rows} odd rows, {bit_minima} bit minima")


def check_proportions():
    examples = 0
    for a in range(1, 11):
        for j in range(1, 8):
            b = 2 * j * a
            for h in range(4, 9):
                t = 2**h
                w = (2 * t + a) // (2 * a)
                parts = (t, t - 1, 2 * j * t)
                assert 2 * abs(t - a * w) <= a
                assert 2 * abs(t - 1 - a * w) <= a + 2
                assert 2 * abs(2 * j * t - b * w) <= b
                coefficient = comb(sum(parts), parts[0]) * comb(parts[1] + parts[2], parts[1])
                assert coefficient % 2 == 1
                examples += 1
    print(f"Bounded-error examples with odd coefficient: {examples}")


if __name__ == "__main__":
    check_weights(192)
    check_proportions()
