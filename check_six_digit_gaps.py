"""Certify weights with six base-p digit units and p dividing N-1."""

from argparse import ArgumentParser
from bisect import bisect_right
from collections import Counter
from functools import lru_cache
from itertools import combinations, combinations_with_replacement
from math import comb, gcd, isqrt, lcm, prod

from flint import fmpz

from check_six_bit_rays import normalized_rows
from search_even_endpoint import avoids, prime_powers


@lru_cache(maxsize=65536)
def factors(n):
    return prime_powers(n)


def candidate_pairings(p, top):
    powers = [p**i for i in range(top+1)]
    for u in range(1, top+1):
        first = powers[u]+1
        common = [gcd(first, value+1) for value in powers]
        gaps = [d for d in range(1, top) if common[d] >= 3]
        for index, d in enumerate(gaps):
            for e in gaps[index:]:
                G = gcd(common[d], common[e])
                if G < 3:
                    continue
                sd, se = common[d]//G, common[e]//G
                ceiling = isqrt(G**3*(sd*se)**2//27)
                span = bisect_right(powers, ceiling)-1
                if e > span:
                    continue
                # Bound the span of u,i,i+d,j,j+e at each extension.
                for i in range(max(1, u-span), min(top-d, u+span-d)+1):
                    start = max(1, max(u, i+d)-span)
                    if d == e:
                        start = max(start, i)
                    stop = min(top-e, min(u, i)+span-e)
                    for j in range(start, stop+1):
                        positions = (u, i, i+d, j, j+e)
                        if p <= 5 and max(positions.count(z) for z in positions) >= p:
                            continue
                        row = (first, powers[i]*(powers[d]+1), powers[j]*(powers[e]+1))
                        n = sum(row)
                        m = (n-1)//powers[min(positions)]
                        if m <= ceiling:
                            s = gcd(m, sd)*gcd(m, se)
                            if 27*(m//s)**2 <= G**3:
                                # The size sieve does not use the spacing theorem.
                                H = max(positions)-min(positions)
                                assert u % d == u % e == 0 and (u//d) % 2 and (u//e) % 2
                                assert H < min(d+e, 3*min(d, e)+gcd(d, e))
                                yield n, row, G, m, s


def assignable(bound, number, n, row):
    # A full power of number must divide bound and its base must be avoided.
    return prod(w for q, w in factors(bound).items()
                if (number//w) % q and avoids(n, row, q))


def complement_by_parts(target, powers):
    # List all single parts, then solve the ordinary integer sum with a lookup.
    modulus = prod(powers)
    residues = [0]
    for w in powers:
        rest = modulus//w
        unit = target*rest*pow(rest, -1, w)
        residues += [(r+unit) % modulus for r in residues]
    parts = {a for r in residues for a in range(r or modulus, target, modulus)}
    return {(a, b, target-a-b) for a in parts for b in parts
            if a <= b <= target-a-b and target-a-b in parts and
            all(sum(x % w == 0 for x in (a, b, target-a-b)) == 2 for w in powers)}


def verify(p, top):
    counts = Counter()
    for n, row, G, m, s in candidate_pairings(p, top):
        counts["size bound"] += 1
        u = assignable(G, n, n, row)
        if u < 3 or 27*(m//s)**2 > u**3:
            continue
        s = assignable(s, m, n, row)
        shift = m//s
        if 27*shift**2 > u**3:
            continue
        counts["after assigned-base carries"] += 1
        shifted = factors(shift)
        if any(2*w+1 > u for w in shifted.values()):
            continue
        counts["CRT searches"] += 1
        powers = tuple(shifted.values())
        rows = list(normalized_rows(u, powers))
        assert {tuple(sorted(y)) for y in rows} == complement_by_parts(u, powers)
        for y in rows:
            counts["CRT rows"] += 1
            v = n//u
            other = tuple(v*z for z in y)
            if all(avoids(n, other, q) for q in shifted):
                if all(avoids(n, other, q) for q in factors(v)):
                    raise ArithmeticError(("candidate requires independent verification", n, row, other))
    print(f"Base {p}, maximum exponent {top}: {dict(counts)}", flush=True)


def brute_pairings(p, top):
    # Enumerate digit strings and partitions independently of the gap inequalities.
    result = set()
    choices = combinations if p == 2 else combinations_with_replacement
    for places in choices(range(1, top+1), 5):
        if max(places.count(i) for i in places) >= p:
            continue
        units = [p**i for i in places]
        n = 1+sum(units)
        m = (n-1)//units[0]
        for t, value in enumerate(units):
            four = units[:t]+units[t+1:]
            for j in (1, 2, 3):
                row = tuple(sorted((1+value, four[0]+four[j], sum(four)-four[0]-four[j])))
                G = gcd(*row)
                cap = gcd(m, lcm(*(gcd(x, y) for x, y in combinations(row, 2))))
                if G >= 3 and 27*(m//cap)**2 <= G**3:
                    result.add((n, row, G, m, cap))
    return result


def cross_check():
    total = 0
    for p in (2, 3, 5, 7):
        generated = [(n, tuple(sorted(row)), G, m, s)
                     for n, row, G, m, s in candidate_pairings(p, 16)]
        assert len(generated) == len(set(generated))
        assert set(generated) == brute_pairings(p, 16)
        total += len(generated)
    # Exact coefficients, including rows which avoid multiple tested primes.
    positives = comparisons = 0
    for n in range(3, 65):
        for a in range(1, n//3+1):
            for b in range(a, (n-a)//2+1):
                row = a, b, n-a-b
                coefficient = comb(n, a)*comb(n-a, b)
                for p in (2, 3, 5, 7):
                    assert avoids(n, row, p) == (coefficient % p != 0)
                    positives += coefficient % p != 0
                    comparisons += 1
    print(f"Matched {total} pairings and {comparisons} exact divisibility checks "
          f"including {positives} prime avoidances.", flush=True)


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("base", type=int, nargs="?")
    parser.add_argument("top", type=int, nargs="?")
    args = parser.parse_args()
    if args.base is None and args.top is None:
        cross_check()
        for base, top in ((2, 159), (3, 63), (5, 39), (7, 39)):
            verify(base, top)
    elif args.top is None or args.top < 1 or not fmpz(args.base).is_prime():
        parser.error("supply a prime base and a positive maximum exponent")
    else:
        verify(args.base, args.top)
