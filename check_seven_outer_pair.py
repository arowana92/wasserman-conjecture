"""Modular-power sieve for the outer-pair gap restrictions."""
from collections import Counter
from fractions import Fraction
from math import gcd, isqrt

from flint import fmpz

BASES = (2, 3, 5, 7, 11, 13, 17, 19, 23)


def prime(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n)+1))


def check_base(p):
    counts = Counter()
    exceptions = []
    limits = (68, 1382) if p == 2 else (26, 899)
    for g in range(limits[0], limits[1]+1):
        if p == 2 and not (prime(g) or g & (g-1) == 0):
            continue
        if p != 2 and g % 3 == 0:
            continue
        # In binary C=1 implies H<3g/2. At odd primes exponent pairs
        # reduce modulo 2g, including negative differences after translation.
        m = (3*g-1)//2 if p == 2 else 2*g-1
        powers = [fmpz(1)]
        for _ in range(max(g, m)):
            powers.append(p*powers[-1])
        odd_f = 1+powers[g]
        while odd_f % 2 == 0:
            odd_f //= 2
        if p == 2:
            low = powers[(4*max(g,101)-3*g)//6]
        else:
            exponent = 4*max(101,3*g)-9*g
            low = powers[max(0, exponent//18-1)]
        for e in range(int(p == 2), m+1):
            # G | J and G | p^g+1 imply (1+p^e)^g = (-1)^(g+f) modulo G.
            w = pow(1+powers[e], g, odd_f)
            counts["powers"] += 1
            for parity in (0, 1):
                cap = odd_f.gcd(w-(-1)**(g+parity))
                if cap <= low:
                    continue
                first = e + int(p == 2)
                for f in range(first+(parity-first)%2, m+1, 2):
                    candidate = cap.gcd(1+powers[e]+powers[f])
                    counts["gcds"] += 1
                    if p == 2:
                        h = max(101,g,f)
                        assert candidate <= powers[(4*h-3*g)//6]
                    elif candidate > low:
                        exceptions.append((g,e,f,int(candidate)))
                        assert 256**6*candidate**18 < 27**6 * fmpz(p)**exponent
    if p == 3:
        assert exceptions == [(58,16,86,12413), (58,30,46,12413), (58,70,100,12413)]
    else:
        assert not exceptions
    expected = {2:70081, 3:1389, 5:224, 7:230}
    assert counts["powers"] == (210225 if p == 2 else 539566)
    assert counts["gcds"] == expected.get(p,0)
    print(f"Base {p}: {counts['powers']} modular powers, {counts['gcds']} direct gcds", flush=True)


def controls():
    assert 3**5 < 2**8
    assert 1383 > 864*Fraction(8,5)
    assert (Fraction(900,6)-1)**3 > 4*900**2
    assert Fraction(27,256) > Fraction(1,27)
    example = (2,12007,33013,55001,65537,65539)
    assert all(prime(e) for e in example) and example[-1]-2 == example[-2]
    for p in BASES:
        for u in range(1, 65):
            for d in range(1, 65):
                g = gcd(u,d)
                common = gcd(1+p**u,1+p**d)
                assert common == (1+p**g if (u//g)%2 and (d//g)%2 else gcd(2,p+1))


if __name__ == "__main__":
    controls()
    for base in BASES:
        check_base(base)
