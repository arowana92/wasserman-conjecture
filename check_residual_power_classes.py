"""Check unrestricted binary/ternary families and power obstructions at sparse layers."""

from functools import cache
from itertools import combinations, product
from math import gcd

from check_endpoint_pair_bound import prime_powers
from check_three_unit_layers import digit_units, free_part, grouped_table


@cache
def full_powers(n):
    return tuple(prime_powers(n))


@cache
def power_residues(J, modulus):
    """A second test, by powering each unit directly."""
    return {pow(t, J, modulus) for t in range(modulus) if gcd(t, modulus) == 1}


def odd_power_unit(value, J, r, direct):
    if direct:
        return value % r in power_residues(J, r)
    return pow(value, (r-1)//gcd(r-1, J), r) == 1


def paired_classes(table, J, direct):
    results = set()
    for p, q, w in table:
        s = w+1
        for U in range(3, w+1):
            if w % U or gcd(U, p*q*s) != 1:
                continue
            # W must equal -1 at each of these odd primes.
            required = {r for r, power in full_powers(s) if 2*power >= U and r != 2}
            required.update(r for r, _ in full_powers(U) if r != 2)
            two_part = U & -U
            for a, b in product(range(J), repeat=2):
                if p == 2:
                    assert s*s >= 4*J
                    if direct:
                        if pow(q, b, 4*J) not in power_residues(J, 4*J):
                            continue
                    elif pow(q, b, 4*J) != 1:
                        continue
                elif not odd_power_unit(q**b, J, p, direct):
                    continue
                if not odd_power_unit(p**a, J, q, direct):
                    continue
                sigma = p**a * q**b
                if not all(odd_power_unit(-sigma, J, r, direct) for r in required):
                    continue
                if two_part > 1:
                    assert sigma % 2 and p > 2
                    modulus = 2*two_part
                    if direct:
                        targets = {t for t in range(modulus)
                                   if (t+1) % two_part == 0 and (s*t+1) % modulus == two_part}
                        if not any(sigma*z % modulus in targets for z in power_residues(J, modulus)):
                            continue
                    else:
                        # The exact 2-part of N fixes its residue modulo 2u.
                        residue = (two_part-1) * pow(s, -1, modulus) % modulus
                        if (residue+1) % two_part or (residue-sigma) % min(modulus, 4*J):
                            continue
                results.add((p, q, w, U, a, b))
    return results


def carry_below(x, y, p, modulus):
    """Legendre layers, including the carry out of the highest inspected digit."""
    place = p
    while place <= modulus:
        if x % place + y % place >= place:
            return True
        place *= p
    return False


def binary_ternary(D):
    # Separate choices of the square root are allowed at different moduli.
    r, = [p for p, _ in full_powers(D)]
    left = set()
    for x in range(1, D):
        if x % r == 0:
            continue
        for ell, _ in full_powers(x):
            if ell <= 3:
                continue
            for a, b in product(range(2), repeat=2):
                sigma = 2**a * 3**b
                if carry_below(x*3**b, (D-x)*3**b, 2, 8):
                    continue
                if carry_below(x*2**a, (D-x)*2**a, 3, 3):
                    continue
                if (odd_power_unit(x*sigma, 2, r, True)
                        and odd_power_unit(-D*sigma, 2, ell, True)):
                    left.add((x, ell, a, b))
    return left


def direct_exponent_periods(D):
    # Independent of the square-character test. Walk both full cyclic orbits
    # modulo the joint congruence modulus and track the two digit cuts.
    def orbit(base, modulus):
        z = 1
        values = []
        while not values or z != 1:
            values.append(z)
            z = base*z % modulus
        return values

    for x in range(1, D):
        for ell, _ in full_powers(x):
            if ell <= 3 or gcd(x, D) > 1:
                continue
            m = D*ell
            for P in orbit(2, 3*m):
                if (x*P) % 3 + ((D-x)*P) % 3 >= 3:
                    continue
                for Q in orbit(3, 8*m):
                    if (x*Q) & ((D-x)*Q) & 7:
                        continue
                    assert x*P*Q % D != 1 or (1+D*P*Q) % ell != 0


def small_static_cofactors(limit):
    # D is a prime power, so D=s. A moving base is at most s-1,
    # since it occurs in a two-unit sum for s-1 with a positive exponent.
    primes = [d for d in range(2, limit) if full_powers(d) == ((d, d),)]
    result = []
    for s in range(3, limit+1):
        if len(full_powers(s)) != 1:
            continue
        w = s-1
        bases = [p for p in primes if p <= w and s % p and len(digit_units(w, p)) == 2]
        for p, q in combinations(bases, 2):
            if free_part(w, (p,q)) >= 3:
                result.append((p, q, w))
    return result


def paired_binary_ternary_cofactors():
    # Both moving bases have order dividing 20 modulo 55.
    allowed = set()
    by_character = set()
    for E in range(1, 1320):
        if gcd(E, 1320) != 1:
            continue
        for a, b in product(range(20), repeat=2):
            if E * 2**a % 3 == 1 and E * 3**b % 8 == 1:
                if E * 2**a * 3**b % 55 == 54:
                    allowed.add(E)
                    break
        if E % 8 in (1, 3):
            a_sign = 1 if E % 3 == 1 else -1
            b_sign = 1 if E % 8 == 1 else -1
            if pow(E, 5, 11) == (-a_sign) % 11 and pow(E, 2, 5) == (a_sign*b_sign) % 5:
                by_character.add(E)
    assert by_character == allowed and min(allowed) == 43 and len(allowed) == 40


def arbitrary_odd_base(D, J):
    r, = [p for p, _ in full_powers(D)]
    for x in range(1, D):
        if x % r == 0 or carry_below(x, D-x, 2, 4*J):
            continue
        for ell, _ in full_powers(x):
            if ell <= 2:
                continue
            for a in range(J):
                if (odd_power_unit(x*2**a, J, r, True)
                        and odd_power_unit(-D*2**a, J, ell, True)):
                    return False
    return True


if __name__ == '__main__':
    for D in (5, 7, 11, 19):
        assert not binary_ternary(D)
        direct_exponent_periods(D)
    print('Unrestricted powers at bases 2 and 3 for cofactors 5, 7, 11 and 19.')
    paired_binary_ternary_cofactors()
    table = grouped_table()
    composite_remainders = {(p, q, w+1) for p, q, w in table if len(full_powers(w+1)) > 1}
    assert composite_remainders == {
        (2, 11, 133), (2, 19, 21), (2, 67, 69), (2, 131, 133),
        (3, 83, 85), (3, 89, 91), (3, 251, 253), (11, 131, 133),
        (13, 47, 2211), (13, 181, 183),
    }
    exceptional = {2: {(2, 17, 18, 3), (2, 17, 18, 9), (2, 71, 72, 3),
                       (2, 71, 72, 9), (2, 191, 192, 3), (5, 29, 30, 3),
                       (5, 149, 150, 3), (13, 181, 182, 7)},
                   4: {(2, 71, 72, 3), (2, 71, 72, 9), (13, 181, 182, 7)},
                   8: {(13, 181, 182, 7)}}
    for J in (2, 4, 8):
        results = paired_classes(table, J, False)
        assert results == paired_classes(table, J, True)
        assert {row[:4] for row in results} == exceptional[J]
        assert all((a % 2, b % 2) == ((0, 1) if p == 5 else (1, 0))
                   for p, q, w, U, a, b in results)
        print(f'J={J}: {len(results)} classes before the companion tests.')
    small = small_static_cofactors(1032)
    assert len(small) == 17 and set(small)-set(table) == {(3,269,270), (3,809,810)}
    assert not paired_classes(small, 8, True)
    remaining = [d for d in range(3, 251, 2) if len(full_powers(d)) == 1]
    previous = set()
    expected = {2: {3,5,7,9}, 4: {17,49,81}, 8: {25,73,89,97,233},
                16: {113,193}, 32: {137}}
    for J in (2, 4, 8, 16, 32):
        good = {D for D in remaining if arbitrary_odd_base(D, J)}
        assert previous <= good and good-previous == expected[J]
        print(f'J={J}: additional universal odd-base cofactors {sorted(good-previous)}.')
        previous = good
