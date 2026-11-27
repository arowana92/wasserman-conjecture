"""Verify prime and residue certificates with Python integer arithmetic."""
from functools import lru_cache
from itertools import combinations_with_replacement
from math import comb, gcd, isqrt, lcm, prod


# (prime divisor of t, multiplicative order of seven)
OBSTRUCTIONS = [
    (19, 3), (43, 6), (117307, 18), (1628413557556843, 54),
    (1357105535093947, 66), (115868130379, 78), (90920746649833, 102),
    (24006520939, 114), (275430251632699963, 114), (1506625093, 126),
    (3426171911854086574033958308360873411, 162),
]

# Distinct prime factors of n-1 for recursive Lucas primality certificates.
FACTORS = {
    3426171911854086574033958308360873411: [2,3,5,7,11,4327,1813629081616156838306060737],
    1813629081616156838306060737: [2,3,13,726614215391088476885441],
    726614215391088476885441: [2,5,7,53,1954488799,3131459623],
    3131459623: [2,3,199,311,937],
    1954488799: [2,3,3121,11597],
    1506625093: [2,3,7,61,98011],
    275430251632699963: [2,3,19,71,809,14020989449],
    14020989449: [2,521,3363961],
    24006520939: [2,3,19,491,22573],
    90920746649833: [2,3,11,17,2879,413923],
    115868130379: [2,3,7,11,13,2756009],
    1357105535093947: [2,3,11,1871,10989954611],
    10989954611: [2,5,193,5694277],
    1628413557556843: [2,3,7,19,37,1063],
}


@lru_cache(None)
def certify_prime(n):
    if n < 10**9:
        assert n >= 2 and all(n % d for d in range(2, isqrt(n)+1))
        return
    fac = FACTORS[n]
    remaining = n-1
    for r in fac:
        certify_prime(r)
        assert remaining % r == 0
        while remaining % r == 0:
            remaining //= r
    assert remaining == 1
    assert any(pow(g, n-1, n) == 1
               and all(gcd(pow(g, (n-1)//r, n)-1, n) == 1 for r in fac)
               for g in range(2, 1000))


def verify_obstruction(r, order):
    certify_prime(r)
    assert pow(7, order, r) == 1
    assert all(pow(7, j, r) != 1 for j in range(1, order))
    period = lcm(6, order)
    modulus = 72*r
    target = 72*(6*pow(72, -1, r) % r)
    even = [pow(7, e, modulus) for e in range(0, period, 2)]
    odd = [7*x % modulus for x in even]
    left = {sum(xs) % modulus for xs in combinations_with_replacement(even, 3)}
    assert all((target-sum(xs)) % modulus not in left
               for xs in combinations_with_replacement(odd, 3))


if __name__ == "__main__":
    for p, order in OBSTRUCTIONS:
        verify_obstruction(p, order)
    short = [sum(7**e for e in exps)
             for exps in combinations_with_replacement(range(6), 6)]
    short = [n for n in short if n % 72 == 0]
    assert len(short) == 33
    assert all(n % 19 != 6 and n % 43 != 6 for n in short)
    h = prod([3,5,7,11,13,17,19,23])
    assert 2*comb(h+5, 6)*10**31 < 2**255
    print(len(OBSTRUCTIONS), "prime obstructions verified")
    print("eight odd primes give exceptional upper density below 1e-31")
