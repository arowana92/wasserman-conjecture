"""Check five digit units at a prime dividing the weight, with bounded full power."""

from argparse import ArgumentParser
from collections import Counter
from itertools import combinations_with_replacement
from math import gcd, isqrt, prod

from flint import fmpz

from search_endpoint import avoids, find_row, smallest_factors


def factors(n):
    result = [(int(p), int(p)**e) for p, e in fmpz(n).factor()]
    assert prod(w for p, w in result) == n and all(fmpz(p).is_prime() for p, w in result)
    return result


def partitions(units):
    """Only the 1+2+2 splits can avoid a predecessor prime."""
    if len(units) != 5:
        return
    seen = set()
    for index, part in enumerate(units):
        others = units[:index] + units[index+1:]
        for j in (1, 2, 3):
            b = others[0]+others[j]
            c = sum(others)-b
            row = part, min(b, c), max(b, c)
            if row not in seen:
                seen.add(row)
                yield row


def check_weight(p, P, units, counts):
    n = sum(units)
    m = n-1
    V = n//P
    if V < 3:
        return
    bounds = {}
    if 27*m*m <= P**3:
        bounds[1] = None  # Also covers all splits containing at least two singletons.
    for a, b, c in partitions(units):
        s = gcd(a-1, b, c)
        if s > 1 and 27*m*m <= P**3*s*s:
            bounds.setdefault(s, []).append((a, b, c))
    if not bounds:
        return
    counts['passed first size bound'] += 1
    shifted = factors(m)
    v_atoms = None
    tests = set()
    for s, rows in bounds.items():
        possible = [(q, w) for q, w in shifted if s % w == 0]
        # All primes available to the first row may be assigned to it.
        if not possible:
            tests.add(1)
            continue
        for row in rows:
            tests.add(prod(w for q, w in possible if avoids(row, q)))
    for S in sorted(tests):
        T = m//S
        required = [(q, w) for q, w in shifted if S % q != 0]
        if 27*T*T > P**3 or len(required) < 3 or any(2*w+1 > P for q, w in required):
            continue
        counts['companion CRT searches'] += 1
        print(f'CRT {p} {P} {n} {S}', flush=True)
        if v_atoms is None:
            v_atoms = factors(V)
        other = find_row(n, V, v_atoms, T, required)
        if other is not None:
            # The S=1 search may have been added as a relaxation even if
            # there is no first row. Report it for exact row-pair review.
            raise ArithmeticError((p, P, n, S, other))
    counts['largest weight passing first size bound'] = max(
        counts['largest weight passing first size bound'], n)


def scan(p, P, counts, prior):
    terms = []
    value = P
    while value < P**3:
        terms.append(value)
        value *= p
    for size in range(3, 6):
        for rest in combinations_with_replacement(terms, size-1):
            units = (P, *rest)
            if p <= size and any(units.count(u) >= p for u in set(units)):
                continue
            n = sum(units)
            if n <= prior or n >= P**3:
                continue
            counts['digit multisets above prior bound'] += 1
            check_weight(p, P, units, counts)


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--limit', type=int, default=1000000)
    parser.add_argument('--prior', type=int, default=10**9,
                        help='skip the previously verified initial interval')
    parser.add_argument('--higher-only', action='store_true',
                        help='restrict to exponents at least two')
    args = parser.parse_args()
    if args.limit < 2 or not 0 <= args.prior <= 10**9:
        parser.error('require limit >= 2 and 0 <= prior <= 1000000000')
    top_prime = isqrt(args.limit) if args.higher_only else args.limit
    spf = smallest_factors(top_prime)
    counts = Counter()
    for p in range(2, top_prime+1):
        if spf[p] != p:
            continue
        P = p*p if args.higher_only else p
        while P <= args.limit:
            counts['full powers'] += 1
            scan(p, P, counts, args.prior)
            P *= p
    print(f'Bound {args.limit}, prior interval through {args.prior}: {dict(counts)}')


if __name__ == '__main__':
    main()
