"""Independent factorial-floor checks for two cuts opposite the same trinomial."""

from argparse import ArgumentParser
from collections import Counter
from itertools import combinations, combinations_with_replacement
from math import comb, gcd, prod

from check_three_factors import powers, valuation


EXAMPLE = 4218421148
SMALL_CUTS = (141136768, 1393201304)
THREE_CUT_WEIGHT = 796691358564


def predecessor_cuts(n, factors):
    m = n - 1
    candidates = [(m//w) * pow(m//w, -1, w) for w in factors.values()]
    return [a for a in candidates if all(valuation(n, (a, n-a), p) == 0 for p in factors)]


def common_scale(n, x, y, factors):
    return prod(w for p, w in factors.items()
                if not valuation(n, (x, n-x), p) and not valuation(n, (y, n-y), p))


def small_scan(bound):
    # Construct powers by trial division, compute each CRT inverse separately,
    # and factor the whole weight for pairs that pass the predecessor tests.
    atoms = []
    for w in range(2, bound+1):
        factors = powers(w)
        if len(factors) == 1:
            atoms.append(next(iter(factors.items())))
    counts = [0, 0, 0]
    for triple in combinations(atoms, 3):
        factors = dict(triple)
        if len(factors) < 3:
            continue
        counts[0] += 1
        n = 1 + prod(factors.values())
        cuts = predecessor_cuts(n, factors)
        if len(cuts) < 2:
            continue
        counts[1] += 1
        weight = powers(n)
        counts[2] += sum(common_scale(n, x, y, weight) >= 3 for x, y in combinations(cuts, 2))
    return counts


def certified_example():
    n = EXAMPLE
    q = n//4
    row = (q, q, 2*q)
    weight = powers(n)
    predecessor = powers(n-1)
    assert weight == {2: 4, q: q}
    assert predecessor == {269: 269, 3137: 3137, 4999: 4999}
    primes = [2, q, *predecessor]
    assert [valuation(n, row, p) for p in primes] == [50, 0, 2, 2, 1]
    cuts = predecessor_cuts(n, predecessor)
    assert sorted(min(a, n-a) for a in cuts) == list(SMALL_CUTS)
    for a in cuts:
        assert [valuation(n, (a, n-a), p) for p in primes] == [0, 1, 0, 0, 0]

    # Confirm modulo each prime directly using Lucas products of small binomials.
    def lucas(total, part, p):
        answer = 1
        while part:
            total, upper = divmod(total, p)
            part, lower = divmod(part, p)
            if lower > upper:
                return 0
            answer = answer * comb(upper, lower) % p
        return answer

    for p in primes:
        assert bool(lucas(n, q, p)*lucas(n-q, q, p) % p) == (p == q)
        assert all(bool(lucas(n, a, p)) == (p != q) for a in cuts)

    oriented = sorted([*cuts, *(n-a for a in cuts)])
    assert all(x & ~y for x, y in combinations(oriented, 2))
    assert all(valuation(n, (x, y-x, n-y), 2) > 0 for x, y in combinations(oriented, 2))
    assert all(sum(row) != n for row in combinations_with_replacement(oriented, 3))
    return n, row, oriented


def minimum_certificate():
    n = EXAMPLE
    q, r = n//4, n-55
    predecessor = powers(n-1)
    assert powers(r) == {r: r} and 2*r > n and n-r <= min(predecessor.values())
    m = n-1
    residues = [(m//w)*pow(m//w, -1, w) for w in predecessor.values()]
    assert sum(residues) == 2*m+1
    witnesses = {
        (2, 269): (5380, 10760),
        (2, 3137): (12548, 50192),
        (2, 4999): (459908, 559888),
        (269, 3137): (1, 843853),
        (269, 4999): (1, 1344731),
        (3137, 4999): (1, 15681863),
    }
    assert set(witnesses) == set(combinations((2, *predecessor), 2))
    for pair, (a, b) in witnesses.items():
        row = (a, b, n-a-b)
        assert min(row) > 0 and all(not valuation(n, row, p) for p in pair)
    # Three rows for which the five-prime list supplies only its external prime.
    rows = [(q,q,2*q), (1,843853,n-843854), (459908,559888,n-1019796)]
    assert all(sum(row) == n and min(row) > 0 for row in rows)
    assert all(any(valuation(n, row, p) == 0 for row in rows) for p in (2,q,*predecessor))
    assert all(valuation(n, row, r) == 1 for row in rows)
    print('Minimum five-prime certificate:', sorted((q,r,*predecessor)))


def check_candidates(filename):
    seen = set()
    weights = set()
    blockers = Counter()
    row_count = 0
    successes = []
    with open(filename) as records:
        for line in records:
            n, u, x, y, a, b, c = map(int, line.split())
            signature = (n, min(x, y), max(x, y))
            assert signature not in seen and x != y
            seen.add(signature)
            weights.add(n)
            predecessor = powers(n-1)
            assert sorted(predecessor.values()) == [a,b,c] and c <= 20000
            cuts = predecessor_cuts(n, predecessor)
            assert x in cuts and y in cuts
            weight = powers(n)
            assert u == common_scale(n, x, y, weight) and u >= 3
            v = n//u
            required = [p for p in weight if v % p == 0]
            rows = [(v*i, v*j, v*(u-i-j))
                    for i in range(1, u//3+1) for j in range(i, (u-i)//2+1)]
            row_count += len(rows)
            blocker = next((p for p in required
                            if all(valuation(n, row, p) > 0 for row in rows)), None)
            if blocker:
                blockers[blocker] += 1
            else:
                successes.extend((n, row, sorted((min(x,n-x), min(y,n-y)))) for row in rows
                                 if all(valuation(n, row, p) == 0 for p in required))
    n, row, _ = certified_example()
    assert (len(seen), len(weights), row_count) == (75, 73, 54435)
    assert blockers == {2: 73, 3: 1}
    assert successes == [(n, row, list(SMALL_CUTS))]
    print('Candidate pairs, weights, normalized rows:', len(seen), len(weights), row_count)
    print('Single-prime exclusions:', dict(blockers))
    print('Successful pair and companion:', successes[0])


def three_singletons():
    n = THREE_CUT_WEIGHT
    predecessor = powers(n-1)
    row = tuple(predecessor_cuts(n, predecessor))
    assert len(row) == 3 and sum(row) == n
    assert gcd(*row) % 31 == 0
    for p in (31, *predecessor):
        assert all(valuation(n, (a,n-a), p) == 0 for a in row)
    assert [valuation(n, row, p) for p in (31, 5521, 11657, 12379)] == [2, 0, 0, 1]


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('candidates', nargs='?', help='output of check_three_factors 20000 --cut-pairs')
    args = parser.parse_args()
    print('Example:', certified_example())
    three_singletons()
    minimum_certificate()
    counts = small_scan(500)
    assert counts == [235224, 2525, 0]
    print('Independent bound-500 counts:', *counts)
    if args.candidates:
        check_candidates(args.candidates)
