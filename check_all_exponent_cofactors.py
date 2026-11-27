"""Build the cofactor and the prime progressions in the uniform delay theorems."""

from fractions import Fraction

from flint import fmpz

from check_cofactor_index_depth import check_247, valuation
from check_prescribed_common_prime import crt
from prime_certificates import primes_through


TARGETS = {2: Fraction(-1, 247), 3: Fraction(-16, 121), 13: Fraction(43), 19: Fraction(127)}
AUXILIARY = {247: 251, 1: 257, 2: 263, 3: 269, 4: 271, 5: 277, 6: 281}


def residue(x, modulus):
    return x.numerator*pow(x.denominator, -1, modulus) % modulus


def least_prime(a, m, lower):
    candidate = a + max(0, (lower-a+m-1)//m)*m
    while not fmpz(candidate).is_prime():
        candidate += m
    return candidate


def prime_cofactor(e):
    factors = {int(p): int(v) for p, v in fmpz(e).factor()}
    equations = []
    for r, target in TARGETS.items():
        modulus = r**(1+int(r == 2)+factors.get(r, 0))
        equations.append((residue(target, modulus), modulus))
    for r in sorted((set(factors)-set(TARGETS)) | {251}):
        equations.append((-1 if r in (31, 251) else 1, r))
    a, m = crt(equations)
    assert m <= 251*4*3*13*19*e
    s = least_prime(a, m, 2)
    assert s > 247
    for r in factors.keys()-TARGETS.keys():
        assert (1+247*s) % r
    return s


def trillion_certificate():
    bound = 10**12
    equations = []
    for r, target in TARGETS.items():
        largest = 1
        while largest*r <= bound:
            largest *= r
        modulus = r**(1+int(r == 2))*largest
        if r == 2:
            two_part = modulus
            modulus *= 2
            target = Fraction(two_part-1, 247)
        equations.append((residue(target, modulus), modulus))
    s, step = crt(equations)
    for index in range(2000000):
        n = 1+247*s
        last = n//(3*two_part)
        assert n == 3*two_part*last and last % 2
        if fmpz(s).is_prime() and fmpz(last).is_prime():
            print(f"Uniform cofactor: s={s}, last weight prime={last}, CRT step index={index}")
            return s, last
        s += step
    raise ArithmeticError("The reported certificate was not found.")


def independent_base_classes():
    cutoff, depth = 47, 8
    equations = []
    for r in primes_through(cutoff):
        if r == 2:
            m = 2**(2*depth)
            target = Fraction(-3*2**depth-1, 247)
        elif r == 3:
            old = 3**depth
            m = 3*old
            low = residue(TARGETS[3], old)
            target = next(Fraction(x) for x in range(low, m, old)
                          if sum(y % m for y in (114*x, 117*x, 16*x+1)) >= m)
        elif r in TARGETS:
            m, target = r**depth, TARGETS[r]
        else:
            m, target = r, Fraction(-1 if r == 31 else 1)
        equations.append((residue(target, m), m))
    for c, p in AUXILIARY.items():
        equations.append(((p-1)*pow(c, -1, p*p) % (p*p), p*p))
    a, m = crt(equations)
    s = least_prime(a, m, 248)
    q = least_prime(1, m, 247*s+2)
    print(f"Independent prime classes: modulus={m}\ns={s}\nq={q}")
    return s, q


if __name__ == "__main__":
    check_247()
    samples = (1, 2, 3, 6, 31, 251, 361, 2756, 2**12*3**5*13**2*19*31*251)
    for e in samples:
        s = prime_cofactor(e)
        for r, target in TARGETS.items():
            m = r**(1+int(r == 2)+valuation(e, r))
            assert residue(target/s, m) == 1
    print(f"Verified prime multipliers at {len(samples)} individual indices.")
    trillion_certificate()
    independent_base_classes()
