"""Periodic certificates that P*t+1 divides its binomial cut at (t^-1 mod P)*t."""

from math import lcm, prod

P = 1009
BASE = 1013
START = 1000
STOP = 11000


def prime_divisors(n):
    factors = []
    divisor = 2
    while divisor * divisor <= n:
        if n % divisor == 0:
            factors.append(divisor)
            while n % divisor == 0:
                n //= divisor
        divisor += 1
    if n > 1:
        factors.append(n)
    return factors


def certificate_depth(exponent, a, prime):
    """Return a local modulus exponent sufficient for all lifts of t."""
    ordinary = weight_valuation = 0
    modulus = 1
    for depth in range(1, 257):
        modulus *= prime
        t = pow(BASE, exponent, modulus)
        total = (P * t + 1) % modulus
        part = a * t % modulus
        ordinary += part > total
        if total:
            if ordinary >= weight_valuation:
                return depth
        else:
            weight_valuation += 1
    raise AssertionError((exponent, a, prime))


def build_certificate():
    factors = [[]] + [prime_divisors(a) for a in range(1, P)]
    depths = {}
    nontrivial = 0
    deepest = None
    inverse = pow(BASE, -1, P)
    a = pow(inverse, START, P)
    for exponent in range(START, STOP):
        for prime in factors[a]:
            depth = certificate_depth(exponent, a, prime)
            depths[prime] = max(depths.get(prime, 0), depth)
            nontrivial += depth > 1
            if deepest is None or depth > deepest[0]:
                deepest = (depth, exponent, prime, a)
        a = a * inverse % P
    modulus = P * prod(p**h for p, h in depths.items())
    period = lcm(P - 1, *((p - 1) * p**(h - 1) for p, h in depths.items()))
    assert pow(BASE, period, modulus) == 1
    assert nontrivial == 7958 and deepest == (37, 5050, 2, 722)
    assert len(depths) == 69 and (modulus.bit_length(), period.bit_length()) == (848, 410)
    return modulus, period


if __name__ == "__main__":
    modulus, period = build_certificate()
    print(f"10,000 exponents, 7,958 nontrivial local checks\nM = {modulus}\nL = {period}")
