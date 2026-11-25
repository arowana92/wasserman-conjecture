"""Verify the lower bounds on certificate sizes in the 3..5000 table."""

from itertools import combinations
from math import comb, gcd, isqrt
from pathlib import Path

# Each row is encoded by its two smallest parts. The first row in each
# sample also needs to be (1,1,n-2) for the core-prime argument below.
SAMPLES = {
    699: [
        (1, 1), (1, 2), (1, 42), (1, 18), (3, 48),
        (3, 3), (3, 279), (3, 174), (84, 246), (3, 36), (18, 324),
        (1, 349), (233, 233),
    ],
    717: [
        (1, 1), (1, 4), (1, 64), (1, 68), (1, 72), (1, 12),
        (1, 8), (1, 196), (1, 76), (9, 192), (3, 3), (3, 93),
        (3, 90), (15, 351), (69, 324), (12, 165), (6, 36),
        (1, 179), (239, 239),
    ],
}


def has_four_for_sample(n, primes):
    """A finite witness against every choice of up to four primes."""
    values = []
    for a, b in SAMPLES[n]:
        if not 0 < a <= b <= n-a-b:
            raise ValueError(f"Invalid row at {n}")
        values.append(comb(n, a) * comb(n-a, b))
    assert values[0] == n*(n-1)
    if gcd(*values) != 1:
        return True
    masks = {p: sum(1 << i for i, value in enumerate(values) if value % p)
             for p in primes if p <= n}
    core = [p for p in masks if values[0] % p == 0]
    # Two primes outside the core are both avoided by the first row.
    for triple in combinations(core, 3):
        if any(masks[p] & masks[q] for p, q in combinations(triple, 2)):
            continue
        used = masks[triple[0]] | masks[triple[1]] | masks[triple[2]]
        if any(used & mask == 0 for mask in masks.values()):
            return True
    return False


def digit_sum(n, p):
    total = 0
    while n:
        n, d = divmod(n, p)
        total += d
    return total


def main():
    path = Path(__file__).with_name("certificates_5000.txt")
    lines = path.read_text().splitlines()
    if len(lines) != 4998:
        raise ValueError("Expected all 4,998 weights")
    primes = [p for p in range(2, 5001)
              if all(p % d for d in range(2, isqrt(p)+1))]
    for n, line in enumerate(lines, start=3):
        weight, *certificate = map(int, line.split())
        if n != weight or len(certificate) not in (1, 4, 5):
            raise ValueError(f"Unexpected entry at {n}")
        if len(certificate) == 1:
            continue  # one is the trivial lower bound
        # Any universal prime divides the coefficient of (1,1,n-2).
        if any(n*(n-1) % p == 0 and digit_sum(n, p) < 3 for p in primes if p <= n):
            raise ValueError(f"A single prime would suffice at {n}")
        if len(certificate) == 5 and has_four_for_sample(n, primes):
            raise ValueError(f"Four primes survived the sample at {n}")
    print("Verified all claimed lower bounds for 3..5000.")


if __name__ == "__main__":
    main()
