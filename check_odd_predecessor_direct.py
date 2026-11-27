"""Reconstruct reported CRT rows independently and verify divisibility by factorial floors."""

from argparse import ArgumentParser
from collections import Counter
from math import comb
from pathlib import Path

from check_three_factors import valuation


def check(rows):
    counts = Counter()
    for p, q, k, x, y, z in rows:
        P = 1 << k
        n = 1+P*p*q
        parts = (P*x, P*y, P*z+1)
        assert sum(parts) == n and min(parts) >= 1
        required = tuple((n-1)//w*pow((n-1)//w, -1, w) for w in (q, p, P))
        assert required == parts and valuation(n, parts, 2) == 0
        vp, vq = (valuation(n, parts, f) for f in (p, q))
        assert vp or vq, (p, q, k, parts)
        counts[bool(vp), bool(vq)] += 1
        if n <= 10000:
            coeff = comb(n, parts[0])*comb(n-parts[0], parts[1])
            assert coeff%2 and (coeff%p == 0) == bool(vp) and (coeff%q == 0) == bool(vq)
    return counts


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("rows", type=Path)
    args = parser.parse_args()
    rows = [tuple(map(int, line.split())) for line in args.rows.read_text().splitlines()]
    assert len(rows) == len({row[:3] for row in rows})
    print("factorial-floor counts at (p,q)", check(rows))
