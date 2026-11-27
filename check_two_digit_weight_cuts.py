"""Check two-unit prefixes against exact coefficients and the earlier residue tests."""

from argparse import ArgumentParser
from itertools import combinations_with_replacement
from math import comb, gcd, prod

from check_predecessor_weight_gap import criterion, factorizations
from check_three_factors import powers
from check_weight_gap_assignments import pair_stage


def least_units(value, p):
    units = []
    place = 1
    while value and len(units) < 3:
        value, digit = divmod(value, p)
        units.extend([place] * min(digit, 3-len(units)))
        place *= p
    return units


def tests(n, weight, predecessor):
    if len(predecessor) < 2:
        return "outside", 0
    q, r = sorted(predecessor.values())[-2:]
    d = (n-1) // (q*r)
    units = {p: least_units(n // power, p) for p, power in weight.items()}
    if any(len(u) < 3 for u in units.values()):
        return "universal weight prime", 0
    profiles = {}
    for p, (_, second, third) in units.items():
        if third > d*(1+second):
            profiles[p] = second
    independent = len(set(profiles.values())) >= 2
    single = False
    for p, ratio in profiles.items():
        power = weight[p]
        others = prod(u[1] for q, u in units.items() if q != p)
        s = 1+ratio
        split_excluded = d <= s*(4*s*others+1)
        grouped_excluded = power == 2 or d % (power*s-1) != 0
        single |= split_excluded and grouped_excluded
    if independent and single:
        return "both tests", d
    if independent:
        return "two primes only", d
    if single:
        return "one prime only", d
    return "inconclusive", 0


def exact_check(bound):
    factors = factorizations(2, bound)
    weights = rows = pairs = 0
    for n in range(3, bound+1):
        _, d = tests(n, factors[n-2], factors[n-3])
        if not d:
            continue
        weights += 1
        divisors = set()
        for a in range(1, n//3+1):
            choose = comb(n, a)
            for b in range(a, (n-a)//2+1):
                value = choose * comb(n-a, b)
                divisors.add(gcd(n*(n-1), value))
                rows += 1
        for x, y in combinations_with_replacement(divisors, 2):
            assert gcd(x, y) > 1, (n, x, y)
            pairs += 1
    print(f"Exact coefficients through {bound}: {weights} weights, {rows} rows, {pairs} gcd pairs.")


def sample_interval(start, count):
    from collections import Counter
    factors = factorizations(start-1, start+count-1)
    groups = Counter()
    additional = Counter()
    for i in range(count):
        n = start+i
        h, m = factors[i+1], factors[i]
        tag, d = tests(n, h, m)
        groups[tag] += 1
        if d and criterion(n, h, m)[0] == "inconclusive":
            additional[tag] += 1
    print(f"Interval of {count} from {start}: {dict(sorted(groups.items()))}")
    print(f"Previously inconclusive inequality tests: {dict(sorted(additional.items()))}")


def example():
    n = 1255603200006
    h, m = powers(n), powers(n-1)
    d = (n-1) // (17*191841589)
    assert d == 385 and prod(h.values()) == n and prod(m.values()) == n-1
    assert n//2 % 2048 == 3 and n//3 % 2187 == 4
    assert 2048 > 3*d and 2187 > 4*d
    assert tests(n, h, m) == ("two primes only", d)
    assert criterion(n, h, m)[0] == "inconclusive" and pair_stage(n, h, m) == 1
    print(f"Two-prime example {n}: cofactor {d}, earlier residue assignments exist.")


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--exact", type=int, default=360)
    parser.add_argument("--start", type=int, default=1001234567890)
    parser.add_argument("--count", type=int, default=100000)
    args = parser.parse_args()
    if args.exact < 3 or args.start < 3 or args.count < 1:
        parser.error("need exact and start >= 3 and count >= 1")
    exact_check(args.exact)
    example()
    sample_interval(args.start, args.count)
