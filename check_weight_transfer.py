"""Finish relaxed endpoint assignments from check_endpoint_range.cpp --transfer."""
import argparse
from collections import Counter
from math import comb, gcd
import sys

from check_range_arithmetic import factor, possible, subsets
from search_endpoint import find_row


def valuation(n, row, p):
    result = 0
    w = p
    while w <= n:
        result += n // w - sum(x // w for x in row)
        w *= p
    return result


def unshifted_free(row, p):
    for j in range(3):
        slopes = list(row)
        slopes[j] -= 1
        if min(slopes) > 0 and valuation(sum(slopes), slopes, p) == 0:
            return True
    return False


def positive_partitions(n):
    for a in range(1, n // 3 + 1):
        for b in range(a, (n - a) // 2 + 1):
            yield a, b, n - a - b


def finish(lines):
    counts = Counter()
    binary_residues = Counter()
    for line in lines:
        n, u, s, x, y, z, missed = map(int, line.split())
        heavy = (x, y, z)
        v, t = n // u, (n - 1) // s
        assert u*v == n and s*t == n-1 and gcd(u, v) == gcd(s, t) == 1
        assert s >= v and sum(heavy) == n and all(a > 0 and a % u == 0 for a in heavy)
        assert all(valuation(n, heavy, p) == 0 for p, _ in factor(s))
        bad = [p for p, _ in factor(u) if valuation(n, heavy, p)]
        assert bad == ([missed] if missed else [])
        counts["heavy rows"] += 1
        counts["largest companion sum"] = max(counts["largest companion sum"], u)
        pred_primes = [p for p, _ in factor(t)]
        scale_primes = [p for p, _ in factor(v)]
        for parts in positive_partitions(u):
            light = tuple(a*v for a in parts)
            counts["companion rows"] += 1
            if any(valuation(n, light, p) for p in pred_primes):
                continue
            light_bad = [p for p in scale_primes if valuation(n, light, p)]
            if len(bad) + len(light_bad) > 1:
                continue
            counts["pairs before shift test"] += 1
            if missed == 2 or light_bad == [2]:
                opposite = light if missed == 2 else heavy
                binary_residues[tuple(sorted(a % 8 for a in opposite))] += 1
            if missed and not unshifted_free(light, missed):
                continue
            if light_bad and not unshifted_free(heavy, light_bad[0]):
                continue
            print("Unresolved specialization:", n, heavy, light, bad + light_bad)
            counts["unresolved"] += 1
    for residues, count in sorted(binary_residues.items()):
        print("Opposite row modulo eight:", residues, "pairs:", count)
    return counts


def independent_window(first, last):
    counts = Counter()
    for n in range(first, last + 1):
        predecessor = factor(n-1)
        for group, s in subsets(predecessor):
            if len(group) < 3:
                continue
            rest = [pair for pair in predecessor if pair not in group]
            for scale, u in subsets(factor(n)):
                v = n // u
                if s < v or not possible(group, s, v) or not possible(rest, (n-1)//s, u):
                    continue
                counts["assignments"] += 1
                for omitted in scale:
                    tested = [pair for pair in scale if pair != omitted]
                    row = find_row(n, u, tested, s, group)
                    if row:
                        counts["passing omissions"] += 1
                        print("Independent search:", n, u, s, row, omitted[0])
    return counts


def controls():
    for n in range(3, 75):
        for row in positive_partitions(n):
            value = comb(n, row[0]) * comb(n-row[0], row[1])
            for p, _ in factor(n*(n-1)):
                assert (valuation(n, row, p) == 0) == (value % p != 0)
                direct = False
                for j in range(3):
                    parts = list(row)
                    parts[j] -= 1
                    if min(parts) > 0:
                        smaller = comb(n-1, parts[0]) * comb(n-1-parts[0], parts[1])
                        direct |= smaller % p != 0
                assert unshifted_free(row, p) == direct
    assert independent_window(3, 20000) == {"assignments": 6084}
    assert independent_window(8284703, 8284703)["passing omissions"] == 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--controls", action="store_true")
    args = parser.parse_args()
    if args.controls:
        controls()
        print("Exact coefficient and independent assignment checks passed.")
    else:
        result = finish(sys.stdin)
        for key in ("heavy rows", "largest companion sum", "companion rows",
                    "pairs before shift test", "unresolved"):
            print(key, result[key])
        if not result["heavy rows"] or result["unresolved"]:
            raise SystemExit(1)
