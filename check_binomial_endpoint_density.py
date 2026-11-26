"""Check the binomial endpoint count, its constants, and finite examples."""

from fractions import Fraction
from math import comb, gcd

from search_endpoint import atoms, smallest_factors


def identities_and_exponents():
    from sympy import expand, symbols
    s, j, d, e, v, c, f = symbols('S J d e V c f')
    difference = v*(d*c + e*f) - j*s*d*e - 1
    factors = (j*s*d - v*f) * (j*s*e - v*c)
    assert expand(factors - (v**2*c*f - j*s) + j*s*difference) == 0
    assert expand(j*s*d*(d*c + e*f) - f*(j*s*d*e + 1) - (j*s*d**2*c-f)) == 0

    first, second = Fraction(21893, 10000), Fraction(2201, 1000)
    assert (first-1)**4 > 2 and (second-1)**6 > 3
    assert first**820 * second**2385 * Fraction(4, 3)**795 < 2**3972
    assert Fraction(1025, 1024)**1000 < 4
    small = -Fraction(1, 3) + Fraction(5, 3)*Fraction(159, 200)
    assert small == Fraction(119, 120) < Fraction(993, 1000)
    assert small + Fraction(8, 3000) < Fraction(997, 1000)


def table(limit):
    hits = {}
    for n in range(3, limit+1):
        groups = {}
        coefficient = 1
        for a in range(1, n//2 + 1):
            coefficient = coefficient * (n-a+1) // a
            groups.setdefault(gcd(coefficient, n*(n-1)), a)
        for g, a in groups.items():
            b = next((b for h, b in groups.items() if gcd(g, h) == 1), None)
            if b is not None:
                assert gcd(n*(n-1), gcd(comb(n, a), comb(n, b))) == 1
                hits[n] = (a, b)
                break
    return hits


def valuation(k, p):
    v = 0
    while k % p == 0:
        v += 1
        k //= p
    return v


def check_assignment(n, a, b, spf):
    ca, cb = comb(n, a), comb(n, b)
    if valuation(ca, 2) > valuation(cb, 2):
        a, b, ca, cb = b, a, cb, ca
    D, J = gcd(n, gcd(ca, cb)), gcd(n-1, gcd(ca, cb))
    u = v = s = d = e = 1
    for p, full in atoms(n, spf):
        left, right = valuation(ca, p), valuation(cb, p)
        retained = full // gcd(full, p**min(left, right))
        if left <= right:
            u *= retained
        else:
            v *= retained
    for p, full in atoms(n-1, spf):
        left, right = valuation(ca, p), valuation(cb, p)
        retained = full // gcd(full, p**min(left, right))
        if left <= right:
            s *= retained
        elif b % retained == 0:
            d *= retained
        else:
            e *= retained
    assert D*u*v == n and J*s*d*e == n-1
    assert a % u == 0 and a*(a-1) % s == 0
    c, r1 = divmod(b, v*d)
    f, r2 = divmod(n-b, v*e)
    assert r1 == r2 == 0 and c > 0 and f > 0
    assert 4*c*f*d*e <= (D*u)**2
    f1, f2 = J*s*d-v*f, J*s*e-v*c
    assert min(f1, f2) > 0 and f1*f2 == v*v*c*f-J*s
    dividend = J*s*d*d*c-f
    assert dividend > 0 and dividend % (d*c+e*f) == 0


if __name__ == '__main__':
    identities_and_exponents()
    hits = table(5000)
    assert sorted(hits) == [46, 370, 476, 477, 609, 974, 989, 2509, 2988, 4324]
    print('Endpoint-coprime binomial witnesses:', hits)
    spf = smallest_factors(5000)
    for n, (a, b) in hits.items():
        check_assignment(n, a, b, spf)
    checks = 0
    for n in range(4, 65):
        for a in range(2, n//2+1):
            for b in range(a, n//2+1):
                check_assignment(n, a, b, spf)
                checks += 1
    print(f'{checks:,} additional partial-valuation assignments passed.')
