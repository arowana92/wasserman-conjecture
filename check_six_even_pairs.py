"""Check repeated binary blocks remaining at six ones and N = 2 modulo 4."""

from argparse import ArgumentParser
from collections import Counter
from itertools import combinations
from math import prod

from check_sparse_companions import factor, eligible_entries
from check_three_factors import valuation


def rows(limit):
    for gap in range(1, limit):
        common = 1 + (1 << gap)
        factors = factor(common)
        for j, k in combinations(range(2, min((gap + 1) // 2, limit - gap)), 2):
            low = 2 + (1 << j) + (1 << k)
            n = common * low
            if 27 * (n - 1)**2 > (2 * common)**3:
                continue
            row = (2 * common, (1 << j) * common, (1 << k) * common)
            u = 2 * prod(w for p, w in factors if low % p and valuation(n, row, p) == 0)
            yield n, row, u


def inspect(limit):
    counts = Counter()
    max_n = 0
    for n, row, u in rows(limit):
        counts['size candidates'] += 1
        assert n.bit_length() <= limit and n % 4 == 2 and valuation(n, row, 2) == 0
        m = n - 1
        if 27 * m * m > u**3:
            continue
        counts['binary row passes size'] += 1
        factors = factor(m)
        if len(factors) < 3 or any(2 * w + 1 > u for _, w in factors):
            continue
        counts['companion CRT systems'] += 1
        max_n = max(max_n, n)
        entries = eligible_entries(n, u, m, factors)
        counts['eligible entries'] += len(entries)
        available = set(entries)
        v = n // u
        for i, a in enumerate(entries):
            for b in entries[i:]:
                c = u - a - b
                if c < b:
                    break
                if c in available:
                    counts['companion triples'] += 1
                    companion = (v * a, v * b, v * c)
                    if all(valuation(n, companion, p) == 0 for p, _ in factors):
                        counts['pass predecessor primes'] += 1
                        if all(valuation(n, companion, p) == 0 for p, _ in factor(v)):
                            raise ArithmeticError((n, row, companion))
    print(limit, dict(counts), 'largest CRT weight', max_n, flush=True)


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('bits', type=int)
    inspect(parser.parse_args().bits)
