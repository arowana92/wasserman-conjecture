"""Independent checks by truncated factorials and direct coefficients."""

from itertools import product
from math import comb, gcd, isqrt, prod

# The verifier uses only the published integers, independently of the generator.
M = 949114553352665995304772175265323965976306160619276984225736539442980999084247742182085764716250905331344992369546515322392183726795472683449657033381251525700344586463220426272624942303193356000041212490463243647373036299737054535422910136320000000000000
L = 1727183973649511851607921419045962918997380239037489529657506994194747518497839483191669351506372578344769159168000000000000


def local_moduli(n):
    local = {}
    for p in range(2, 1010):
        power = 1
        while n % p == 0:
            n //= p
            power *= p
        if power > 1:
            local[p] = power
    assert n == 1
    return local


def floor_sum(n, parts, p, limit):
    total = 0
    m = p
    while m <= limit:
        total += n // m - sum(x // m for x in parts)
        m *= p
    return total


def check_block():
    local = local_moduli(M)
    assert pow(1013, L, M) == 1 and local[1009] == 1009
    prime = 166 * M + 1
    divisors = set(local) | {2, 83}
    assert all(all(p % i for i in range(2, isqrt(p) + 1)) for p in divisors | {1013})
    assert pow(181, prime - 1, prime) == 1
    assert all(gcd(pow(181, (prime - 1) // p, prime) - 1, prime) == 1 for p in divisors)
    t = pow(1013, 1000, M)
    needed = 0
    for _ in range(10000):
        # A full set of residues suffices, irrespective of the chosen lift.
        a = next(x for x in range(1, 1009) if (x * t - 1) % 1009 == 0)
        n = 1009 * t + 1
        b = a * t
        remaining = a
        for p, limit in local.items():
            count = 0
            while remaining % p == 0:
                remaining //= p
                count += 1
            if count and n % p == 0:
                needed += 1
                # b*C(n,b) = n*C(n-1,b-1), with v_p(b)=v_p(a).
                assert floor_sum(n - 1, (b - 1, n - b), p, limit) >= count
        assert remaining == 1
        t = t * 1013 % M
    assert needed == 7958
    return needed


def check_digit_counts():
    columns = layers = 0
    for prime in (2, 3, 5, 7):
        for r in (1, 2):
            passing = sum(sum(xs[:-1]) <= xs[-1]
                          for xs in product(range(prime), repeat=r + 1))
            assert passing == comb(prime + r, r + 1)
            columns += prime ** (r + 1)
            # Compare the digit witness to the actual multinomial floor.
            nmax = 36
            for xs in product(range(1, nmax), repeat=r):
                for n in range(sum(xs) + 1, nmax + 1):
                    parts = xs + (n - sum(xs),)
                    lower = 1
                    while lower <= n:
                        m = prime * lower
                        exceeds = sum(x % m // lower for x in xs) > n % m // lower
                        cost = n // m - sum(x // m for x in parts)
                        assert cost >= int(exceeds)
                        lower = m
                        layers += 1
    return columns, layers


def check_direct():
    cases = failures = 0
    for p in range(2, 33):
        for t in range(1, 40):
            if gcd(t, p) != 1:
                continue
            n = p * t + 1
            a = pow(t, -1, p)
            b = a * t
            coeff = comb(n, b)
            assert coeff * b == n * comb(n - 1, b - 1)
            if coeff % n:
                failures += 1
                assert gcd(a, n) > 1
            cases += 1
    assert failures > 0
    n = 177996
    parts = (22020, 13580, 142396)
    assert sum(parts) == n and gcd(*parts) == 4
    assert [floor_sum(n, (x, n - x), 2, n) for x in parts] == [1, 0, 3]
    assert floor_sum(n, parts, 2, n) == 9
    assert prod(p**e for p, e in ((2, 2), (3, 1), (7, 1), (13, 1), (163, 1))) == n
    assert all(floor_sum(n, parts, p, n) >= e
               for p, e in ((2, 2), (3, 1), (7, 1), (13, 1), (163, 1)))
    return cases, failures


def check_power_representatives():
    cases = 0
    for b in range(2, 80):
        for p in range(2, 160):
            if gcd(b, p) != 1:
                continue
            powers = {1}
            high = b
            while high < p:
                powers.add(high)
                high *= b
            value = b % p
            fits = True
            while value != 1:
                if value not in powers:
                    fits = False
                    break
                value = value * b % p
            assert fits == ((high - 1) % p == 0)
            cases += 1
    p = (10**19 - 1) // 9
    factors = (2, 3, 5, 7, 11, 13, 19, 37, 52579, 333667)
    assert all(all(n % i for i in range(2, isqrt(n) + 1)) for n in factors)
    assert 3 * prod(factors) == p - 1
    assert pow(3, p - 1, p) == 1
    assert all(gcd(pow(3, (p - 1) // q, p) - 1, p) == 1 for q in factors)
    order = (p - 1) // 6
    assert pow(2, order, p) == 1
    assert all(pow(2, order // q, p) != 1 for q in factors if order % q == 0)
    assert pow(10, 19, p) == 1 and 10 % p != 1
    stride = order // gcd(order, 19)
    assert stride == 9746588693957115
    mersenne = 2**31 - 1
    partner = mersenne + 2**20
    assert partner == 2148532223
    assert all(all(n % i for i in range(2, isqrt(n) + 1)) for n in (mersenne, partner))
    orbit = {10**j for j in range(19)}
    for u in range(1, 25):
        for j in (-7, -1, 0, 1, 6):
            v = u + j * stride
            if v <= 0:
                continue
            t = pow(2, u, p) * pow(5, v, p) % p
            assert pow(t, -1, p) in orbit
    return cases, stride


if __name__ == "__main__":
    print("Nontrivial shifted factorial checks:", check_block())
    print("Digit columns and Legendre layers:", check_digit_counts())
    print("Exact binomials and full-weight failures:", check_direct())
    print("Power representative checks and repunit stride:", check_power_representatives())
