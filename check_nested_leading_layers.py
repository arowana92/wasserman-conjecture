"""Certify leading base-prime columns for nested predecessor powers."""

from fractions import Fraction
from math import gcd, lcm

from check_three_factors import powers

# (D, U, x, t, k): (base, denominator C for E/W <= 1/C).
CASES = {
    (5, 3, 3, 1, 1): (2, 4),
    (7, 3, 3, 3, 1): (2, 1),
    (7, 3, 3, 2, 2): (2, 1),
    (7, 5, 5, 1, 1): (2, 6),
    (8, 3, 3, 4, 1): (2, 1),
    (8, 5, 5, 2, 1): (2, 1),
    (9, 4, 4, 4, 1): (2, 1),
    (9, 5, 5, 3, 1): (3, 2),
    (9, 5, 5, 2, 2): (2, 4),
    (10, 3, 3, 6, 1): (3, 2),
    (10, 3, 3, 5, 2): (3, 1),
    (10, 3, 3, 3, 4): (3, 1),
    (10, 3, 6, 3, 1): (3, 2),
    (10, 3, 6, 2, 2): (2, 2),
    (10, 7, 7, 2, 1): (5, 9),
}


def floor_failure(w, slopes, layer):
    return sum(int(a*w/layer) for a in slopes[:3]) < int(slopes[3]*w/layer)


def leading_cover(x, t, k, r, C):
    D = x+t+k
    coeff = (Fraction(x), Fraction(t), Fraction(C*k+1, C), Fraction(D))
    cuts = {Fraction(1), Fraction(r)}
    for j in range(4):
        for a in coeff:
            for h in range(1, int(r*a/r**j)+1):
                point = h*r**j/a
                if 1 < point < r:
                    cuts.add(point)
    endpoints = sorted(cuts)
    for lo, hi in zip(endpoints, endpoints[1:]):
        assert any(floor_failure(lo, coeff, r**j) for j in range(4)), (lo, hi)

    # Independent integral grid, including every discontinuity as a left endpoint.
    M = lcm(x, t, C*k+1, D)
    for z in range(M, r*M):
        assert any(x*z//w + t*z//w + (C*k+1)*z//(C*w) < D*z//w
                   for w in [M*r**j for j in range(4)])
    return len(cuts)-1, (r-1)*M


def check():
    extras = {(7, 4, 4, 2, 1), (9, 7, 7, 1, 1), (10, 7, 7, 1, 2)}
    candidates = set()
    counts = [0, 0]
    for D in (5, 7, 8, 9, 10):
        b = 5 if D == 10 else D
        for x in range(3, D-1):
            for u in range(3, x+1):
                if x % u or gcd(x, b) != 1 or gcd(D, u) != 1:
                    continue
                for k in range(1, D-x):
                    if gcd(k, b*u) == 1:
                        candidates.add((D, u, x, D-x-k, k))
    assert candidates == set(CASES) | extras
    for (D, u, x, t, k), (r, C) in CASES.items():
        b = 5 if D == 10 else D
        assert r in powers(b*u) or (r == 2 and u < 7)
        assert all((k*P+1) % (b*u) for P in range(2, C))
        for i, n in enumerate(leading_cover(x, t, k, r, C)):
            counts[i] += n
    assert all((119*q % 25 + (16*q+1) % 25 + 35*q % 25) >= 25
               for q in range(4, 25, 5))
    print(f"{len(candidates)} parameter cases, {len(CASES)} leading covers, "
          f"{counts[0]} rational intervals and {counts[1]} grid points")
    print("The remaining P=17 branch carries at the modulus 25 in all five lifts.")


if __name__ == '__main__':
    check()
