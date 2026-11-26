"""Independent exhaustive binomial-product audit through the first mixed example."""

from argparse import ArgumentParser
from collections import defaultdict
from math import comb, gcd


def all_pairs(n):
    endpoint = n * (n - 1)
    by_gcd = defaultdict(list)
    value = 1
    for a in range(1, n//2 + 1):
        value = value * (n+1-a) // a
        by_gcd[gcd(endpoint, value)].append(a)

    found = set()
    for g, cuts in by_gcd.items():
        # Each individual part of a compatible row gives a compatible binomial.
        positions = [b for k, indices in by_gcd.items() if gcd(g, k) == 1 for b in indices]
        pool = set(positions) | {n-b for b in positions}
        positions.sort()
        for i, b in enumerate(positions):
            if 3*b > n:
                break
            for c in positions[i:]:
                d = n-b-c
                if d < c:
                    break
                if d in pool and gcd(g, comb(n-b, c)) == 1:
                    found.update((a, b, c, d) for a in cuts)
    return found


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--limit', type=int, default=5752)
    args = parser.parse_args()
    for n in range(3, args.limit+1):
        found = all_pairs(n)
        if found:
            print(n, sorted(found), flush=True)
