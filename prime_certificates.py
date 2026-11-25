"""Find finite prime certificates for positive trinomial coefficients."""

from argparse import ArgumentParser
from functools import lru_cache
from itertools import combinations
from math import isqrt


def primes_through(limit):
    sieve = bytearray([1]) * (limit + 1)
    sieve[:2] = b'\0\0'
    for p in range(2, isqrt(limit) + 1):
        if sieve[p]:
            sieve[p*p::p] = bytes((limit - p*p) // p + 1)
    return [p for p in range(2, limit + 1) if sieve[p]]


def free_count(n, p):
    """Number of ordered positive carry-free triples, by inclusion-exclusion."""
    three = two = 1
    while n:
        n, d = divmod(n, p)
        three *= (d+1)*(d+2)//2
        two *= d+1
    return three - 3*two + 3


def free_rows(n, p):
    """Allocate digits from most to least significant, keeping a <= b <= c."""
    digits = []
    while n:
        n, d = divmod(n, p)
        digits.append(d)

    def visit(i, a, b, c):
        if i < 0:
            if a:
                yield a, b, c
            return
        d = digits[i]
        for x in range(d+1):
            if a == b and 2*x > d:
                break
            for y in range(x if a == b else 0, d-x+1):
                z = d-x-y
                if b == c and y > z:
                    break
                yield from visit(i-1, a*p+x, b*p+y, c*p+z)

    yield from visit(len(digits)-1, 0, 0, 0)


def free_in_base(row, p):
    a, b, c = row
    while c:
        a, x = divmod(a, p)
        b, y = divmod(b, p)
        c, z = divmod(c, p)
        if x+y+z >= p:
            return False
    return True


def find_cover(n, primes):
    """Return up to three positive triples collectively avoiding every prime."""
    @lru_cache(None)
    def witness(group):
        if not group:
            return (1, 1, n-2)
        p, *others = group  # primes are ordered by increasing free_count
        return next((row for row in free_rows(n, p)
                     if all(free_in_base(row, q) for q in others)), None)

    def partition(i, groups):
        if i == len(primes):
            return [witness(group) for group in groups]
        p = primes[i]
        for j, group in enumerate(groups):
            enlarged = group + (p,)
            if witness(enlarged):
                result = partition(i+1, groups[:j] + [enlarged] + groups[j+1:])
                if result is not None:
                    return result
            if not group:  # other empty groups are equivalent
                break
        return None

    return partition(0, [(), (), ()])


def four_prime_certificate(n, ranked):
    """With no universal prime, search all four mutually disjoint avoidance sets."""
    # (1,1,n-2) avoids every prime outside the prime divisors of n*(n-1).
    # At most one such outside prime can occur in a disjoint family.
    primes = [p for _, p in ranked]
    core = [p for p in primes if n*(n-1) % p == 0]
    counts = {p: count for count, p in ranked}

    @lru_cache(None)
    def disjoint(p, q):
        if counts[p] > counts[q]:
            p, q = q, p
        return not any(free_in_base(row, q) for row in free_rows(n, p))

    for triple in combinations(core, 3):
        if all(disjoint(p, q) for p, q in combinations(triple, 2)):
            for q in primes:
                if q not in triple and all(disjoint(p, q) for p in triple):
                    return tuple(sorted((*triple, q)))
    return ()


def find_certificate(n, candidates):
    ranked = sorted((free_count(n, p), p) for p in candidates if p <= n)
    if not ranked[0][0]:
        return (ranked[0][1],)
    chosen = set()
    while True:
        ordered = [p for _, p in ranked if p in chosen]
        rows = find_cover(n, ordered)
        if rows is None:
            # Delete redundant primes while preserving the certificate.
            for p in reversed(ordered):
                trial = [q for q in ordered if q != p and q in chosen]
                if find_cover(n, trial) is None:
                    chosen.remove(p)
            if len(chosen) > 4:
                smaller = four_prime_certificate(n, ranked)
                if smaller:
                    return smaller
            return tuple(sorted(chosen))
        for _, p in ranked:
            if all(not free_in_base(row, p) for row in rows):
                chosen.add(p)
                break
        else:
            raise ArithmeticError(f"Candidate gcd-one rows at {n}: {rows}")


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("first", type=int)
    parser.add_argument("last", type=int)
    args = parser.parse_args()
    if not 3 <= args.first <= args.last:
        parser.error("require 3 <= first <= last")
    primes = primes_through(args.last)
    for n in range(args.first, args.last+1):
        witness = find_certificate(n, primes)
        print(n, *witness)


if __name__ == "__main__":
    main()
