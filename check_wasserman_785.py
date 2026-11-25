"""Check Wasserman's trinomial conjecture at total weight 785 using Python 3."""

from collections import Counter
from itertools import combinations_with_replacement
from math import comb, gcd


def main():
    weight = 785
    modulus = 2 * 3 * 5 * 7
    divisors = Counter()
    for a in range(1, weight // 3 + 1):
        for b in range(a, (weight - a) // 2 + 1):
            # The remaining part is c = weight - a - b, so 1 <= a <= b <= c.
            coefficient = comb(weight, a) * comb(weight - a, b)
            divisors[gcd(coefficient, modulus)] += 1

    # Every coefficient occurs among the sorted triples by symmetry.
    for x, y, z in combinations_with_replacement(divisors, 3):
        if gcd(gcd(x, y), z) == 1:
            raise RuntimeError(f"No common prime established for divisors {(x, y, z)}")

    print(f"Checked {sum(divisors.values())} unordered positive triples.")
    print(f"Frequencies by gcd with {modulus}: {dict(sorted(divisors.items()))}")
    print("Every three coefficients share a prime from {2, 3, 5, 7}.")


if __name__ == "__main__":
    main()
