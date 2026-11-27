"""Check both prime certificates, actual weight layers, and the exponent guarantees."""

from fractions import Fraction
from math import gcd

from check_scaled_row_blocks import certify_prime
from prime_certificates import primes_through


STATIC = 23128931932329398448536814804548764654513479119879993
LAST = 865967223229951653409533606573487133219537
S = 1306232213153877917593933367823536683066159885737114363997573716868895663450681
Q = 322880541893073397389196821002188355954113056549133075089295162936641374385930241


def first_layer(parts, p, depth):
    total = sum(parts)
    w = 1
    for k in range(1, depth+1):
        w *= p
        if total//w != sum(x//w for x in parts):
            return k
    return None


def trillion():
    for n in (STATIC, LAST):
        certify_prime(n)
    assert 1+247*STATIC == 3*2**41*LAST
    missed = []
    for r, target in [(2, Fraction(-1,247)), (3, Fraction(-16,121)), (13, 43), (19, 127)]:
        difference = Fraction(target, STATIC)-1
        assert difference.denominator % r
        numerator = abs(difference.numerator)
        power = 1
        while numerator % (r*power) == 0:
            power *= r
        missed.append(power//2 if r == 2 else power)
    missed.append(LAST*(LAST-1))
    assert min(missed) > 10**12
    density_loss = sum(Fraction(1, m) for m in missed)
    assert density_loss < Fraction(2024, 10**15)
    print("Sufficient coverage omits at most", float(density_loss), "of all exponents.")


def rectangle():
    for n in (S, Q):
        certify_prime(n)
    for c, p in [(1,257),(2,263),(3,269),(4,271),(5,277),(6,281),(247,251)]:
        assert Q % (p*p) == 1 and (1+c*S) % (p*p) == p
    for p in primes_through(47):
        depth = 16 if p == 2 else 9 if p == 3 else 8 if p in (13, 19) else 1
        assert Q % (p**depth) == 1
    assert Q > 247*S+1
    for b in (1,2,6,19,52,247):
        w = S*Q**b
        n = 1+247*w
        assert n % 9 in (3,6)
        rows = ((1,1,n-2), (114*w,117*w,16*w+1), (n//3,)*3)
        assert all(sum(row) == n for row in rows)
        for p in primes_through(47):
            assert any(first_layer(row, p, 8) is None for row in rows)
        assert [first_layer(row, 3, 9) for row in rows] == [1,9,1]
        assert [first_layer(row, 2, 9) for row in rows] == [1,1,9]
        assert n % 65536 == 64768 and n//3 % 65536 == 65280
        for row in rows[1:]:
            partial = sum(n//(1 << j) - sum(x//(1 << j) for x in row) for j in range(1, 17))
            assert partial >= 8
        assert [first_layer(row, 251, 2) for row in rows] == [1,1,2]
    print("Six exact weights passed the floor checks, including the", n.bit_length(), "bit weight at exponent 247.")


def coset_counts():
    # Directly test the small-representative property at accessible field sizes.
    # For primes below the asymptotic threshold, count the gaps without a size claim.
    cases = covered = 0
    for p in primes_through(4000):
        if p <= 247:
            continue
        for e in (1,2,3,4,6,8,12,24):
            g = gcd(e, p-1)
            image_labels = {pow(z, (p-1)//g, p) for z in range(1, (p-1)//247+1)}
            cases += 1
            covered += len(image_labels) == g
    # At these larger primes the character bound proves coverage for each test.
    for p in (10000019, 10000079, 20000003):
        certify_prime(p)
        for e in (2,3,6):
            g = gcd(e, p-1)
            h = (p-1)//494
            assert h*h > (g-1)**2*p
            assert len({pow(z, (p-1)//g, p) for z in range(2, 2*h+1)}) == g
    print(f"Small representatives cover all cosets in {covered}/{cases} control cases, and all nine large-field tests.")


if __name__ == "__main__":
    trillion()
    rectangle()
    coset_counts()
