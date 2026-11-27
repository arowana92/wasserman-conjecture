"""Check the entropy margins, lattice count and small binary pattern tables."""

from argparse import ArgumentParser
from collections import Counter
from math import gcd


def entropy_below(a, b, v, w):
    """Check H_2(a/b) < v/w by integer powers."""
    return b ** (b * w) < 2 ** (b * v) * (a**a * (b - a) ** (b - a)) ** w


def margins():
    # gamma=85/109, kappa=1+gamma.
    assert entropy_below(79, 342, 85, 109)
    assert entropy_below(79, 114, 97, 109)
    # The positive root of x^3=x^2+1 is below 47/32 < 2^(5/9).
    assert 47**3 > 47**2 * 32 + 32**3 and 47**9 < 2**5 * 32**9
    assert 194 * 9 > 16 * 109 and 194 * 50 < 89 * 109


def lattice_checks():
    tests = 0
    for k in range(2, 100):
        for q in range(2, 36):
            # Full sup-norm shortest vector, including vectors on the axes.
            shortest = k
            for u in range(1, k + 1):
                r = u * q % k
                shortest = min(shortest, max(u, min(r, k - r)))
            for x in (3, 9, 19):
                count = Counter((q * i + j) % k for i in range(x) for j in range(x))
                # Euclidean diameter <= sqrt(2) X, and ||v||_2 <= sqrt(2)L.
                # The two line and spacing bounds give these safe constants.
                rhs = (k + 2 * x * shortest) * (shortest + 2 * x)
                assert max(count.values()) * k * shortest <= rhs
                tests += 1
    return tests


def free(a, b, c):
    return not ((a & b) | (a & c) | (b & c))


def pattern_counts(limit):
    counts = Counter()
    for d in range(1, limit + 1):
        uppers = [
            (a, b, d - a - b)
            for a in range(1, d // 3 + 1)
            for b in range(a, (d - a) // 2 + 1)
            if free(a, b, d - a - b)
        ]
        for c in range(1, limit + 1):
            kind = 0
            for x in uppers:
                for slot in (-1, 0, 1, 2):
                    shift = int(slot >= 0)
                    if any((c - shift) * a % d for a in x):
                        continue
                    y = [(c - shift) * a // d + int(i == slot) for i, a in enumerate(x)]
                    if not free(*y):
                        continue
                    if shift:
                        g = gcd(*x)
                        other = [v // g for i, v in enumerate(x) if i != slot]
                        if gcd(*other) < 3:
                            continue
                    kind |= 1 if slot < 0 else 2
            counts[kind] += 1
    return [counts[j] for j in range(4)]


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--limit", type=int, default=63)
    args = parser.parse_args()
    margins()
    print("Lattice bounds verified on", lattice_checks(), "boxes and lattices")
    print(args.limit, *pattern_counts(args.limit))
