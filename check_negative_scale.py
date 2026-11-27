"""Check a 447-digit prime and the negative-scale carry calculation."""
from math import comb, gcd, prod

from prime_certificates import primes_through


def prime_powers(d):
    result = {}
    for p in primes_through((d+2)//2):
        w = p
        while w*p <= d+1:
            w *= p
        result[p] = w
    return result


def trace(index, initial, modulus):
    """Compute z^index + z^(-index) from z + z^(-1) = initial."""
    a, b = 2, initial
    for digit in bin(index)[2:]:
        middle = (a*b-initial) % modulus
        if digit == "0":
            a, b = (a*a-2) % modulus, middle
        else:
            a, b = middle, (b*b-2) % modulus
    return a


def prime_certificate(initial):
    powers = prime_powers(2026)
    m = 12*prod(powers.values())
    q = m-1
    assert 2 in powers and 3 in powers and len(str(q)) == 447
    if trace(m, initial, q) != 2:
        return False
    return all(gcd(trace(m//p, initial, q)-2, q) == 1 for p in powers)


def carry_controls():
    cases = exact = 0
    for d in range(7, 101):
        for p, w in prime_powers(d).items():
            for multiplier in (1, 3):
                scale = multiplier*w-1
                layer = p*w
                for a in range(p, d, p):
                    for b in range(a, d-a, p):
                        assert a*scale % layer + b*scale % layer == 2*layer-a-b >= layer
                        cases += 1
                        if (a+b)*scale <= 700:
                            assert comb((a+b)*scale, a*scale) % p == 0
                            exact += 1
    print("Forced carries:", cases, "exact binomial checks:", exact)


if __name__ == "__main__":
    carry_controls()
    assert prime_certificate(45001160)
    print("Certified the 447-digit prime at cofactor bound 2026.")
