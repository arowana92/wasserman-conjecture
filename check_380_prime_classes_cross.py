"""Independent prefix, orbit and lifted-root checks for the cofactor 380 sieve."""
from fractions import Fraction
from math import gcd

from check_cofactor_986_powers_cross import extendible_prefix
from prime_certificates import primes_through


def verify():
    slopes = (130, 247, 3)
    for p, k, digit, target in [(5, 4, 3, 458), (13, 2, 4, 56)]:
        assert [x for x in range(digit, p**k, p) if extendible_prefix(slopes, x, p, k)] == [target]
    assert not extendible_prefix(slopes, 234, 19, 2)
    for digit in (1, 3):
        assert sum((a*digit+b) % 4 for a, b in zip(slopes, (0, 0, 1))) == 5

    # Build the full orbit by repeated multiplication; invert generator
    # powers to identify the exceptional exponent for each candidate base.
    m = 2006875
    orbit = []
    w = 1
    while not orbit or w != 1:
        orbit.append(w)
        w = w * 1471708 % m
    length = len(orbit)
    assert length == 58500
    logs = {q: pow(j, -1, length) for j, q in enumerate(orbit) if gcd(j, length) == 1}
    assert len(logs) == 14400 and [logs[p] for p in (1733, 1973, 5717)] == [51541, 52433, 42971]
    for q, e in logs.items():
        assert pow(q, e, m) == 1471708
    # Enumerate the entire binary power orbit modulo 32 in each odd class
    # of exponents modulo 4. Check low bits by addition of the three parts.
    for e in (1, 3):
        good = 0
        for q in range(1, 32, 2):
            for k in range(e, e+32, 4):
                v = (380 * pow(q, k, 32) + 1) * pow(13, -1, 128) % 128
                parts = [a*v for a in (1, 4, 8)]
                if sum(parts) == (parts[0] | parts[1] | parts[2]):
                    good += 1
                    break
        assert good == 5

    # Lift the equation x^PERIOD=1 by its derivative, then test the carry
    # using a difference of factorial floors on an integral companion row.
    mass = Fraction(3, 1976)
    count = 0
    for p in primes_through(100000):
        if (380*length) % p == 0:
            continue
        root = pow(-380, -1, p)
        if pow(root, length, p) != 1:
            continue
        derivative = length * pow(root, length-1, p) % p
        layer = p
        for h in range(2, 6):
            next_layer = layer*p
            defect = (pow(root, length, next_layer)-1) // layer
            root += layer * (-defect*pow(derivative, -1, p) % p)
            layer = next_layer
            assert pow(root, length, layer) == 1
            # Choose an ordinary positive lift so the companion scale is integral.
            w = next(x for x in range(root, root+13*layer, layer) if (380*x+1) % 13 == 0)
            v = (380*w+1)//13
            if 13*v//layer > sum(a*v//layer for a in (1, 4, 8)):
                mass *= 1 - Fraction(1, layer-layer//p)
                count += 1
                break
    assert count == 53 and Fraction(998536633, 10**9) < 1-mass < Fraction(998536636, 10**9)
    print("Independent orbit and lifted-root checks:", len(logs), "classes,", count, "companion primes.")


if __name__ == "__main__":
    verify()
