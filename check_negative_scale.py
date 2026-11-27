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


def check_quadratic_order(initial):
    primes = list(prime_powers(2026))
    m = 12*prod(prime_powers(2026).values())
    q = m-1

    # Pairs are coefficients of 1,z in the quotient by z²-initial*z+1.
    def multiply(x, y):
        a, b = x
        c, d = y
        return (a*c-b*d) % q, (a*d+b*c+initial*b*d) % q

    def power(x, exponent):
        result = (1, 0)
        while exponent:
            if exponent % 2:
                result = multiply(result, x)
            x = multiply(x, x)
            exponent //= 2
        return result

    # On a group of primes the supplied element is z^(m/product(group)).
    def check_group(element, group):
        if len(group) == 1:
            a, b = element
            assert gcd(q, a-1, b) == 1
            return
        middle = len(group)//2
        left, right = group[:middle], group[middle:]
        check_group(power(element, prod(right)), left)
        check_group(power(element, prod(left)), right)

    z = (0, 1)
    assert power(z, m) == (1, 0)
    check_group(power(z, m//prod(primes)), primes)


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
    check_quadratic_order(45001160)
    print("Certified the 447-digit prime at cofactor bound 2026.")
