"""Verify cyclic cut comparisons, finite group filters and binary class counts."""

from functools import lru_cache
from math import gcd

from check_geometric_exponent_pairs import equation_count, multiplier_count, prime_power
from check_independent_fermat_exponents import prime_factors
from check_prime_bases_at_geometric_exponents import prime_by_order


Q = 2**61 - 1


def digit_sum(x, p):
    total = 0
    while x:
        x, digit = divmod(x, p)
        total += digit
    return total


@lru_cache(maxsize=128)
def powers(a, m):
    assert gcd(a, m) == 1
    s = {1}
    t = a % m
    while t not in s:
        s.add(t)
        t = t*a % m
    return s


def passes_group(x, d, ell):
    m = d*ell
    r = pow(x, -1, d)
    target = r + d*(((-pow(d, -1, ell)-r)*pow(d, -1, ell)) % ell)
    assert x*target % d == 1 and d*target % ell == ell-1
    h = powers(2, m)
    k = powers(Q*Q % m, m)
    match = any(target*t % m in h for t in k)
    # Independent target reconstruction and intersection in the other direction.
    other = ell - pow(d, -1, ell)
    other += ell*((r-other)*pow(ell, -1, d) % d)
    assert match == any(other*t % m in k for t in h)
    return match


def tables():
    totals = [0, 0, 0, 0]
    empty = []
    for d in range(3, 251):
        if d % 2 and prime_power(d):
            cuts = {}
            for x in range(3, d):
                if (x & (d-x)) or gcd(x, d) != 1:
                    continue
                primes = [p for p in prime_factors(x) if p != 2]
                if not primes:
                    continue
                key = min(x, d-x)
                cuts[key] = cuts.get(key, False) or any(passes_group(x, d, p) for p in primes)
            if not any(cuts.values()):
                empty.append(d)
            for x, fixed_q in cuts.items():
                count = multiplier_count(x, d-x, 60)
                assert count == equation_count(x, d-x, 60)
                totals[0] += 1
                totals[1] += count
                if fixed_q:
                    totals[2] += 1
                    totals[3] += count
        if d == 100:
            assert totals == [165, 11376825559964, 88, 5222790089452]
            print(f"D <= 100: {totals}")
    assert totals == [760, 35935947104515, 492, 21688449846682]
    assert empty == [3, 5, 7, 9, 17, 31, 49, 97]
    print(f"D <= 250: {totals}; empty at {empty}")


def rotations():
    checked = 0
    for p in (2, 3, 5, 7):
        for d in range(3, 512):
            if gcd(d, p) != 1:
                continue
            orbit = powers(p, d)
            h = len(orbit)
            c = (p**h-1)//d
            sums = [digit_sum(i, p) for i in range(d+1)]
            weight = digit_sum(c, p)
            for x in orbit:
                if sums[x] + sums[d-x] != sums[d]:
                    continue
                assert digit_sum(c*x, p) == weight
                assert h*(p-1)-weight <= (sums[d]-sums[x])*weight
                if (sums[d]-1)*weight < h*(p-1):
                    assert sums[x] == 1
                checked += 1
    for h, q in ((9, 2**17-1), (31, Q)):
        d = 73 if h == 9 else 2**31-1
        prime_by_order(d)
        prime_by_order(q)
        assert (4*q*q-1) % d == 0
        assert (q+1) % (1 << d.bit_length()) == 0
        c = (2**h-1)//d
        assert (digit_sum(d, 2)-1)*digit_sum(c, 2) < h
    print(f"{checked} submasks in power orbits checked")


if __name__ == "__main__":
    rotations()
    tables()
