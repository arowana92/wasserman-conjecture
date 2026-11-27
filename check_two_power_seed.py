"""Finite reference-weight test for two prime-power factors with large exponents."""

from argparse import ArgumentParser
from functools import lru_cache
from itertools import combinations
from math import comb, gcd

from check_three_factors import powers
from prime_certificates import free_in_base, free_rows


def free(row, p):
    return free_in_base(tuple(sorted(row)), p)


def grouped(d):
    """Return a projected binomial cut and companion, if they pass the relaxation."""
    n = d+1
    weight, static = powers(n), powers(d)
    for u in range(3, d):
        if n % u or gcd(u, n//u) > 1:
            continue
        v = n//u
        for x in range(u, d, u):
            cut = (0, x, n-x)
            if any(not free(cut, p) for p in weight if u % p == 0):
                continue
            needed = [p for p in static if not free(cut, p)]
            if any(2*static[p]+1 > u for p in needed):
                continue
            for a in range(1, u//3+1):
                for b in range(a, (u-a)//2+1):
                    parts = (a, b, u-a-b)
                    row = tuple(v*f for f in parts)
                    if any(not free(row, p) for p in needed):
                        continue
                    missing = [p for p in weight if v % p == 0 and not free(row, p)]
                    if len(missing) <= 2 and all(free((0, x, d-x), p) for p in missing):
                        yield x, u, parts, tuple(missing)


def prime_part(value, p):
    w = 1
    while value % p == 0:
        w *= p
        value //= p
    return w


@lru_cache(maxsize=512)
def pair_costs(total, p):
    """Least allowed p-part of the multiplier for each primitive positive pair."""
    triples = list(free_rows(total, p))
    triples.extend((0, i, total-i) for i in range(1, total//2+1)
                   if free((0, i, total-i), p))
    costs = {}
    for row in triples:
        for i, j in combinations(row, 2):
            if i == 0:
                continue
            g = gcd(i, j)
            pair = (i//g, j//g)
            value = prime_part(g, p)
            costs[pair] = min(costs.get(pair, value), value)
    return costs


def minimum_sum(tables, d):
    if not tables:
        return 2
    tables.sort(key=len)
    best = d+1
    for pair, value in tables[0].items():
        cost = sum(pair) * value
        for table in tables[1:]:
            cost *= table.get(pair, d+1)
            if cost >= best:
                break
        best = min(cost, best)
    return best


def split_costs(d):
    """Smallest pair-sum products, capped per side, with and without missing primes."""
    n = d+1
    weight = powers(n)
    primes = tuple((*weight, *powers(d)))
    regular = [pair_costs(n//weight.get(p, 1), p) for p in primes]
    shifted = [pair_costs(d, p) if p in weight else None for p in primes]
    full = (1 << len(primes))-1
    sums = []
    for mask in range(full+1):
        slots = [i for i in range(len(primes)) if (mask >> i) & 1]
        ordinary = minimum_sum([regular[i] for i in slots], d)
        relaxed = ordinary
        for j in slots:
            if shifted[j] is not None:
                relaxed = min(relaxed, minimum_sum(
                    [shifted[i] if i == j else regular[i] for i in slots], d))
        sums.append((ordinary, relaxed))
    return tuple(min(sums[m][j] * sums[full^m][j] for m in range(full+1))
                 for j in range(2))


def literal_split(d):
    """Independent direct scan of four multipliers using exact binomial products."""
    n = d+1
    weight = powers(n)
    primes = tuple((*weight, *powers(d)))
    full = (1 << len(primes))-1

    @lru_cache(None)
    def possible(total, p, a, b, ppart):
        for k in range(1, total//(a+b)+1):
            if prime_part(k, p) == ppart:
                value = comb(total, k*a) * comb(total-k*a, k*b)
                if value % p:
                    return True
        return False

    pairs = []
    for x in range(1, d//4+1):
        for y in range(x, (d-1)//2-x+1):
            g = gcd(x, y)
            ordinary = shifted = 0
            for i, p in enumerate(primes):
                args = (p, x//g, y//g, prime_part(g, p))
                if possible(n//weight.get(p, 1), *args):
                    ordinary |= 1 << i
                if p in weight and possible(d, *args):
                    shifted |= 1 << i
            pairs.append((x+y, ordinary, shifted))
    pairs.sort()
    relaxed_pass = False
    for i, (cost, ok, alt) in enumerate(pairs):
        for j in range(i, len(pairs)):
            other, good, shift = pairs[j]
            if cost*other >= d:
                break
            missing = full ^ (ok|good)
            if missing == 0:
                return True, True
            options = [1 << k for k in range(len(primes)) if (missing >> k) & 1]
            if len(options) == 1:
                relaxed_pass |= bool(missing & (alt|shift))
            elif len(options) == 2:
                left, right = options
                relaxed_pass |= bool((alt & left and shift & right)
                                     or (alt & right and shift & left))
    return False, relaxed_pass


def literal_grouped(d):
    n = d+1
    for u in range(3, d):
        v, remainder = divmod(n, u)
        if remainder or gcd(u, v) != 1:
            continue
        for x in range(u, d, u):
            cut = comb(n, x)
            if gcd(u, cut) > 1:
                continue
            shifted = comb(d, x)
            for a in range(1, u//3+1):
                for b in range(a, (u-a)//2+1):
                    coeff = comb(n, v*a) * comb(n-v*a, v*b)
                    missing = gcd(v, coeff)
                    if (gcd(d, gcd(cut, coeff)) == 1 and gcd(shifted, missing) == 1
                            and len(powers(missing)) <= 2):
                        return True
    return False


def check(limit, cross):
    unrestricted = []
    coprime = []
    grouped_weights = []
    for d in range(1, limit+1):
        cost0, cost = split_costs(d)
        if d <= cross:
            assert literal_split(d) == (cost0 < d, cost < d), d
        first = next(grouped(d), None)
        if d <= cross or d == 475:
            assert literal_grouped(d) == bool(first), d
        if first:
            grouped_weights.append(d)
        if first or cost < d:
            unrestricted.append(d)
        if cost0 < d or (first and any(not g[-1] for g in grouped(d))):
            coprime.append(d)
    print(f"Through {limit}: certified for all bases {limit-len(unrestricted)}, "
          f"for bases coprime to D+1 {limit-len(coprime)}.")
    print("Unresolved cofactors for unrestricted bases:", *unrestricted)
    print("Cofactors with a surviving grouped reference pair:", *grouped_weights)
    print(f"Exact four-multiplier cross-check through {min(cross, limit)}.")


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--limit", type=int, default=1500)
    parser.add_argument("--cross", type=int, default=150)
    args = parser.parse_args()
    if args.limit < 1 or args.cross < 0:
        parser.error("need a positive limit and a nonnegative cross-check bound")
    check(args.limit, args.cross)
