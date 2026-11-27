"""Decide the intermediate gaps using exact prime powers and the heavy CRT row."""

from argparse import ArgumentParser
from collections import Counter
from functools import lru_cache
from math import prod

from flint import fmpz
from search_endpoint import find_row, splits


@lru_cache(maxsize=20000)
def factors(value):
    result = [(int(p), int(p**v)) for p, v in fmpz(value).factor()]
    assert prod(w for _, w in result) == value
    assert all(fmpz(p).is_prime() for p, _ in result)
    return result


def check_weight(n, counts):
    weight, predecessor = factors(n), factors(n - 1)
    if len(weight) < 2 or len(predecessor) < 3:
        return
    for u, first, second in splits(weight):
        v = n // u
        if min(u, v) < 3:
            continue
        for s, left, right in splits(predecessor):
            t = (n - 1) // s
            if s < v or len(left) < 3:
                continue
            if 2 * max(w for _, w in left) >= v or 2 * max((w for _, w in right), default=0) >= u:
                continue
            if 27 * t * t > u**3 or 27 * s * s > v**3:
                continue
            counts["assignments"] += 1
            if find_row(n, u, first, s, left):
                counts["heavy rows"] += 1
                if find_row(n, v, second, t, right):
                    raise ArithmeticError(("unresolved endpoint pair", n, u, s))


def scan(first, last, cap):
    counts = Counter()
    largest = 0
    for d in range(first, last + 1):
        for c in range(1, cap + 1):
            q = 2
            while True:
                n = d * q + c
                if q > c and (d * c) ** 5 < 27 * (n - 1) ** 2:
                    break
                if n > 10**9:
                    counts["weights"] += 1
                    largest = max(largest, n)
                    check_weight(n, counts)
                q *= 2
    return counts, largest


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--first", type=int, default=1)
    parser.add_argument("--last", type=int, default=255)
    parser.add_argument("--cap", type=int, default=255)
    args = parser.parse_args()
    print(scan(args.first, args.last, args.cap))
