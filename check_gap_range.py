"""Verify the non-endpoint case for 5,001 <= N <= 1,000,000 using exact integers."""

from collections import Counter
from functools import lru_cache
from itertools import combinations
from math import gcd, lcm, prod
from pathlib import Path

from prime_certificates import primes_through


@lru_cache(None)
def layer_constants(a, b):
    # For w > max(a,b), the two low residues already sum to a+b >= 11.
    constants = [60]*6
    for w in range(7, max(a, b) + 1):
        for h in range(a % w + b % w, 6):
            constants[h] = lcm(constants[h], w)
    return tuple(constants)


@lru_cache(None)
def groups_for_gap(g):
    groups = {}
    for small_sum in range(11, g + 1):
        for a in range(1, small_sum // 2 + 1):
            b = small_sum - a
            key = (layer_constants(a, b), 1 if a <= 2 else 64)
            groups.setdefault(key, []).append((a, b))
    return groups


def factorial_valuation(n, p):
    total = 0
    while n:
        n //= p
        total += n
    return total


def valuation(n, row, p):
    return factorial_valuation(n, p) - sum(factorial_valuation(x, p) for x in row)


def subnumbers(n, p):
    """Integers with every base-p digit at most the corresponding digit of n."""
    values = [0]
    place = 1
    while n:
        n, digit = divmod(n, p)
        values = [v + j*place for j in range(digit + 1) for v in values]
        place *= p
    return values


def digit_cost(n, p):
    # Counts all ordered carry-free triples, allowing zero parts.
    ways = 1
    while n:
        n, digit = divmod(n, p)
        ways *= (digit + 1)*(digit + 2)//2
    return ways


@lru_cache(None)
def incompatible(n, p, q):
    if digit_cost(n, q) < digit_cost(n, p):
        p, q = q, p
    for a in subnumbers(n, p):
        if not 1 <= a <= n//3:
            continue
        for b in subnumbers(n-a, p):
            if a <= b <= (n-a)//2 and valuation(n, (a, b, n-a-b), q) == 0:
                return False
    return True


def valid_certificate(n, row, witness, primes):
    if len(set(witness)) != len(witness):
        return False
    if any(p not in primes or valuation(n, row, p) == 0 for p in witness):
        return False
    if len(witness) == 1:
        p, = witness
        digit_sum, value = 0, n
        while value:
            value, digit = divmod(value, p)
            digit_sum += digit
        return digit_sum <= 2
    return len(witness) == 3 and all(incompatible(n, p, q) for p, q in combinations(witness, 2))


def main():
    path = Path(__file__).with_name("gap_certificates_1m.txt")
    witnesses = {}
    for line in path.read_text().splitlines():
        n, a, b, *primes = map(int, line.split())
        key = n, a, b
        if key in witnesses:
            raise ArithmeticError(("duplicate certificate", key))
        witnesses[key] = tuple(primes)
    last = 1000000
    primes = set(primes_through(last))
    counts = Counter()
    r = 2
    for n in range(3, last + 1):
        if n in primes:
            r = n
        if n <= 5000:
            continue
        gap = n - r
        if gap < 11:
            counts["gap below 11"] += 1
            continue
        norm = n**18
        for (constants, factor), pairs in groups_for_gap(gap).items():
            lower = prod(((n-h)//gcd(n-h, k))**(6-h) for h, k in enumerate(constants))
            if lower * 3**18 > factor * norm:
                counts["product exclusions"] += len(pairs)
                continue
            for a, b in pairs:
                witness = witnesses.pop((n, a, b), ())
                if not valid_certificate(n, (a, b, n-a-b), witness, primes):
                    raise ArithmeticError(("unproved row", n, a, b, witness))
                counts["prime certificates"] += 1
    if witnesses:
        raise ArithmeticError(("unused certificates", len(witnesses)))
    print(dict(counts))
    print("Excluded non-endpoint counterexamples at all weights from 5001 through 1000000.")


if __name__ == "__main__":
    main()
