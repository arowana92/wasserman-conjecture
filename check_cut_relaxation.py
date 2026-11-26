"""Search for a binomial and a trinomial with no common endpoint prime."""

from argparse import ArgumentParser
from collections import Counter
from math import comb, gcd

from flint import fmpz

from check_central_cuts import binomial_carry
from search_endpoint import atoms, find_row, smallest_factors, splits


def cuts_for_assignment(n, u, weight_atoms, s, predecessor_atoms):
    v = n // u
    residues = [0]
    for _, w in predecessor_atoms:
        rest = s // w
        unit = v * rest * pow(rest, -1, w)
        residues += [(value + unit) % s for value in residues]
    for residue in residues:
        for y in range(residue or s, v // 2 + 1, s):
            a = u * y
            if a >= 2 and all(not binomial_carry(n, a, p) for p, _ in
                              weight_atoms + predecessor_atoms):
                yield a


def cut_witnesses(n, spf):
    found = {}
    weight, predecessor = atoms(n, spf), atoms(n-1, spf)
    if len(weight) < 2 or len(predecessor) < 2:
        return found
    predecessors = list(splits(predecessor))
    for u, cut_weight, row_weight in splits(weight):
        v = n // u
        if u < 3 or v < 2:
            continue
        for s, cut_predecessor, row_predecessor in predecessors:
            t = (n-1) // s
            if 27 * t*t > u**3 or any(2*w >= u for _, w in row_predecessor):
                continue
            cuts = cuts_for_assignment(n, u, cut_weight, s, cut_predecessor)
            first = next(cuts, None)
            if first is None:
                continue
            row = find_row(n, v, row_weight, t, row_predecessor)
            if row is not None:
                found[first] = row
                for a in cuts:
                    found[a] = row
    return found


def direct(n, a, row):
    b, c, d = row
    assert a >= 2 and 2*a <= n and min(row) >= 1 and b+c+d == n
    binomial = fmpz.bin_uiui(n, a)
    trinomial = fmpz.bin_uiui(n, b) * fmpz.bin_uiui(n-b, c)
    assert binomial.gcd(n*(n-1)).gcd(trinomial) == 1


def first_example():
    n, a = 5752, 568
    row = (719, 719, 4314)
    direct(n, a, row)
    trinomial = comb(n, 719) * comb(n-719, 719)
    b = comb(n, a)
    assert gcd(n*(n-1), trinomial) == 46008
    assert gcd(n*(n-1), b) == 719
    counts = Counter()
    for side in (a, n-a):
        # All unordered positive splits of the chosen side.
        split = 1
        for j in range(1, side//2 + 1):
            split = split * (side+1-j) // j
            g = gcd(n*(n-1), gcd(b*split, trinomial))
            counts[tuple(p for p in (2, 3, 71, 719) if g % p == 0)] += 1
    assert counts == {(2, 3): 9, (2, 71): 10, (3, 71): 10, (2, 3, 71): 2847}
    return counts


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--limit', type=int, default=1000000)
    args = parser.parse_args()
    if args.limit < 3:
        parser.error('use a limit of at least 3')
    spf = smallest_factors(args.limit)
    hits = {}
    for n in range(3, args.limit + 1):
        witnesses = cut_witnesses(n, spf)
        if witnesses:
            for a, row in witnesses.items():
                direct(n, a, row)
            hits[n] = sorted(witnesses)
            print(n, hits[n], flush=True)
    print(f'{len(hits)} weights through {args.limit:,}.')
    for n, smaller_cuts in hits.items():
        allowed = set(smaller_cuts) | {n-a for a in smaller_cuts}
        triples = {(a, b, n-a-b) for a in smaller_cuts for b in smaller_cuts
                   if a <= b <= n-a-b and n-a-b in allowed}
        if triples:
            print('All three cuts occur:', n, sorted(triples))
    print('Common-prime sets for the refinements at 5752:', first_example())
