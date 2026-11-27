"""Check entropy margins and exact pairs for separate weight and predecessor divisors."""

from fractions import Fraction as F
from itertools import combinations_with_replacement
from math import comb, gcd

from check_binary_windows import ln_bounds
from prime_certificates import free_rows


def parameters():
    eps = F(1, 10000)
    eta = F(1, 10**8)
    for p, u, alpha in ((2, F(1, 5), F(5245, 10000)),
                        (3, F(29, 100), F(5024, 10000))):
        mean = (p-1)*(F(1, 2) + 2*eps)/3
        logp_low, logp_high = ln_bounds(F(p))
        entropy = (ln_bounds(sum(u**d for d in range(p)))[1]
                   + mean * ln_bounds(1/u)[1]) / logp_low
        assert entropy + 2*alpha/3 < 1 - 2*eta
        assert 2*eps**2 / logp_high > eta
    h_upper = (ln_bounds(F(4))[1] + 3*ln_bounds(F(4, 3))[1])/(4*ln_bounds(F(2))[0])
    assert h_upper < F(8113, 10000)
    print("The binary and ternary entropy bounds and H(1/4) pass rational enclosures.")


def first_pairs(bases, stop):
    first = [None, None]
    eligible = 0
    for n in range(3, stop+1):
        candidates = set().union(*(free_rows(n, p) for p in bases))
        if n <= 90:
            direct = set()
            for a in range(1, n//3+1):
                for b in range(a, (n-a)//2+1):
                    k = comb(n, a) * comb(n-a, b)
                    if any(k % p for p in bases):
                        direct.add((a, b, n-a-b))
            assert direct == candidates
        divisors = {}
        for a, b, c in sorted(candidates):
            if a < 2:
                continue
            k = comb(n, a) * comb(n-a, b)
            assert any(k % p for p in bases)
            eligible += 1
            divisors.setdefault(gcd(k, n*(n-1)), (a, b, c))
        for (g, row), (h, partner) in combinations_with_replacement(divisors.items(), 2):
            for j in (0, 1):
                if first[j] is None and gcd(n-j, g, h) == 1:
                    first[j] = (n, row, partner)
    return first, eligible


if __name__ == "__main__":
    parameters()
    for primes, end, expected in (((2,), 511, (511, 127)),
                                  ((2, 3), 216, (216, 53))):
        examples, count = first_pairs(primes, end)
        assert tuple(x[0] for x in examples) == expected
        print(f"Bases {primes}, {count} interior rows through {end}: {examples}")
