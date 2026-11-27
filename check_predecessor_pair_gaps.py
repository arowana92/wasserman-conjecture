"""Check the square-root window and inverse-residue bounds on predecessor pairs."""

from itertools import combinations
from math import comb, gcd

from check_static_power_transfer import prime_powers
from check_three_factors import valuation


def compositions(s):
    for a in range(s+1):
        for b in range(s-a+1):
            yield a, b, s-a-b


def check_nearby():
    powers = prime_powers(160)
    evaluated = rows = zeros = reconstructed = 0
    for P, Q in combinations(sorted(powers), 2):
        p, q, h = powers[P], powers[Q], Q-P
        if p == q or h > 12:
            continue
        for D in range(1, 13):
            n = D*P*Q + 1
            for a in range(D+1):
                for b in range(D*h+1):
                    if not (a or b):
                        continue
                    norm = a*a + D*b*(b-a*h)
                    scale = gcd(n, a*P+b)
                    assert norm % scale == 0
                    if norm == 0:
                        assert (D, h) in ((1, 2), (4, 1))
                        zeros += 1
                    if P > D*h and scale > 2*Q and norm:
                        quotient = (a*P+b)//scale
                        assert 1 <= quotient <= D and a > 0
                        assert divmod(quotient*scale-b, a) == (P, 0)
                        reconstructed += 1
                    evaluated += 1
            if D > 4 or D*h >= P or D*h*h+1 >= Q:
                continue
            # Exhaust actual parts at the two P-blocks, then use factorial valuations.
            for high in compositions(D):
                for low in compositions(D*h):
                    row = tuple((a*P+b)*P+(i == 2)
                                for i, (a, b) in enumerate(zip(high, low)))
                    if min(row) < 1 or gcd(*row) < 3:
                        continue
                    assert valuation(n, row, p) > 0 or valuation(n, row, q) > 0
                    rows += 1
    print(f"nearby powers: {evaluated} divisor checks, {zeros} zero resultants, "
          f"{rows} actual row pairs of valuation tests, {reconstructed} reconstructions")


def check_linear_window():
    powers = prime_powers(75)
    checked = 0
    for P, Q in combinations(sorted(powers), 2):
        p, q = powers[P], powers[Q]
        if p == q:
            continue
        for D in range(2, 24):
            for base, variable in ((q, P), (p, Q)):
                place = base
                while place <= base*D*variable:
                    if D*((-variable) % place) <= place:
                        assert all(comb(D*variable, x*variable) % base == 0
                                   for x in range(1, D))
                        checked += 1
                        break
                    place *= base
    print(f"Negative-residue obstructions checked against exact binomials: {checked}")


def check_direction(P, p, Q, d, L):
    # This upper addition is necessary even using just one of the two primes.
    occupied = []
    z = Q
    i = 0
    while z:
        z, digit = divmod(z, p)
        if digit:
            occupied.append(i)
        i += 1
    cap = P
    while cap < P*d:
        cap *= p
    assert all(p**(j-i) <= cap for i, j in zip(occupied, occupied[1:]))
    H = p
    while H < Q:
        m, h = divmod(Q, H)
        if d*h < H:
            l, c = divmod(L, H)
            I = m*c-l*h
            J, rem = divmod(P*I-m, Q)
            assert rem == 0 and abs(I) <= d*m*h and abs(J)*Q <= m*(d*P*h+1)
            if m % p:
                assert H < d*P*h+1
            if H == P:
                assert (l*h*h-m*c*h-m*m) % Q == 0
                if Q > d*m*h*h+m*m:
                    assert m*m % h == 0
        H *= p


def check_distinct():
    powers = prime_powers(110)
    probed = accepted = 0
    for P, Q in combinations(sorted(powers), 2):
        p, q = powers[P], powers[Q]
        if p == q:
            continue
        for D in range(2, 12):
            n = 1+D*P*Q
            for x in range(1, D):
                d = D-x
                for L in range(pow(P, -1, Q), d*Q+1, Q):
                    row = (x*P*Q, P*L, n-x*P*Q-P*L)
                    probed += 1
                    if valuation(n, row, p) or valuation(n, row, q):
                        continue
                    check_direction(P, p, Q, d, L)
                    check_direction(Q, q, P, d, row[2]//Q)
                    accepted += 1
    assert accepted
    print(f"distinct positions: {accepted} carry-free rows among {probed} CRT rows")


if __name__ == "__main__":
    check_nearby()
    check_linear_window()
    check_distinct()
