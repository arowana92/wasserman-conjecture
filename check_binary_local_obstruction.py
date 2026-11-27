"""Construct simultaneous low-digit witnesses on the cofactor 12031."""
from math import gcd, lcm

from prime_certificates import primes_through


D = 12031
SLOPES = (681, 9222, 2128)
SEED = 106675
STEP = 22211280
TARGETS = {3: 17140160990152208, 53: 33, 227: 8}
SKIP = (5, 7, 13, 113, 1093, 3511)


def row(q):
    x, y, z = SLOPES
    return x*q, y*q, z*q+1


def valuation(parts, p):
    total = sum(parts)
    w = p
    v = 0
    while w <= total:
        v += total//w - sum(k//w for k in parts)
        w *= p
    return v


def check_seeds():
    x, y, z = SLOPES
    assert D == 53*227 == x+y+z and (x&y | y&z | z&x) == 0
    assert STEP == 16*27*5*7*13*113
    assert SEED % 5876 == 907 and STEP % 5876 == 0
    for r in SKIP:
        assert STEP % (r-1) == 0 and (D*pow(2, SEED, r)+1) % r
    for p, q in TARGETS.items():
        m = p
        while STEP % m == 0:
            m *= p
        assert pow(2, SEED, m) == q % m
        assert pow(2, STEP, m) == 1 and pow(2, STEP, m*p) != 1
        assert valuation(row(q), p) == 0
    assert (D*TARGETS[3]+1) % 9 == 6
    assert [r for r in primes_through(10**6) if r > 2 and pow(2,r-1,r*r) == 1] == [1093,3511]


def lift(a, period, p, target, depth):
    """Impose 2^a = target modulo p^depth while preserving the old progression."""
    modulus = p
    while modulus < p**depth and pow(2, period, modulus*p) == 1:
        modulus *= p
    assert pow(2, a, modulus) == target % modulus and pow(2, period, modulus) == 1
    while modulus < p**depth:
        next_modulus = p*modulus
        value = pow(2, a, next_modulus)
        derivative = (pow(2, period, next_modulus)-1)//modulus
        digit = ((target-value)//modulus)*pow(value*derivative, -1, p) % p
        a += digit*period
        period *= p
        modulus = next_modulus
    return a, period


def witness(primes, depth):
    a, period = SEED, STEP
    for p, q in TARGETS.items():
        a, period = lift(a, period, p, q, depth+1)
    for r in sorted(set(primes) - {2, *TARGETS, *SKIP}):
        if pow(2, r-1, r*r) == 1 or period % r == 0:
            raise ValueError(("prime outside the induction hypotheses", r))
        period = lcm(period, r-1)
        if (D*pow(2,a,r)+1) % r == 0:
            modulus = r**(depth+1)
            target = (3*r-1)*pow(D, -1, modulus) % modulus
            a, period = lift(a, period, r, target, depth+1)
    return a, period


def check_witness(e, primes, depth):
    # Work directly with the actual row residues, including C via division by 3.
    active = []
    for p in primes:
        modulus = p**(depth+1)
        q = pow(2, e, 3*modulus)
        n = D*q+1
        assert n % 3 == 0
        if p > 2 and n % p and (n-1) % p:
            continue  # The full endpoint coefficient is a unit.
        if p != 2 and p not in TARGETS:
            assert n % (p*p) == 3*p % (p*p)
            active.append(p)
        parts = (1, 1, n-2), row(q), (n//3,)*3
        assert any(all(sum(k % (p**i) for k in triple) < p**i
                           for i in range(1, depth+2)) for triple in parts)
    return active


def first_full_predecessor():
    # At this exponent the odd row covers every prescribed lowest layer.
    q, n = 2**119, 1903*2**119+1
    parts = (264*q, 519*q, 1120*q+1)
    assert sum(parts) == n and n % 9 == 3
    assert valuation(parts, 2) == 0 and valuation(parts, 3) == 81
    assert all(sorted(v % w for v in parts) == [0, 0, 1] for w in (q, 11, 173))
    assert all(v % 3 == 0 for v in parts)
    assert sum(v % 9 for v in parts) >= 9
    assert 1903 == 11*173 and all(pow(2,2580,w) == 1 for w in (9,11,173))


if __name__ == "__main__":
    check_seeds()
    first_full_predecessor()
    for bound, depth in [(200, 32), (4096, 12)]:
        primes = list(primes_through(bound))
        a, h = witness(primes, depth)
        for e in (a, a+h):
            active = check_witness(e, primes, depth)
        print(f"Primes through {bound}, depth {depth} beyond the first layer, "
              f"period of {h.bit_length()} bits, additional weight primes {active}.")
