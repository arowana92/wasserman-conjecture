"""Exact finite degree calculations for the seven-unit affine-line theorem."""
from collections import defaultdict
from fractions import Fraction
from itertools import combinations_with_replacement, permutations
from math import gcd

from flint import fmpz_poly
from check_binomial_trinomial_gcd import doubled_offsets
from check_quadratic_digit_positions import phase_tuples
from check_quadrinomial_small_bases import CASES

X = fmpz_poly([0, 1])
P = 3


def block(slopes, powers):
    """Remove unit monomials and clear negative coefficient exponents."""
    start, low = min(slopes), min(powers)
    return sum(P**(a-low) * X**(j-start) for j, a in zip(slopes, powers))


def mixed_quadruples():
    order = CASES[P][0][-1]
    assert order == 12
    offsets = defaultdict(set)
    for heights, count in CASES[P][1].items():
        if count:
            for oriented in set(permutations(heights)):
                delta = tuple(a-oriented[0] for a in oriented)
                for phase in phase_tuples(P, oriented):
                    offsets[phase].add(delta)
    high = set()
    allocations = full = 0
    for n in range(1, 25):
        roots = [k for k in range(order) if k*n % order == order//2]
        for tail in combinations_with_replacement(range(25), 3):
            if gcd(n, *tail) != 1:
                continue
            v = (0, *tail)
            height = max(n, v[-1])
            possible = set()
            for k in roots:
                possible.update(offsets.get(tuple(k*j % order for j in v), ()))
            for delta in possible:
                for parity in (0, 1):
                    doubled = tuple(parity*j+h for j, h in zip(v, delta))
                    if parity*n % 2 or any(a % 2 for a in doubled):
                        continue
                    powers = tuple(a//2 for a in doubled)
                    labels = tuple(zip(v, powers))
                    if any(labels.count(label) >= P for label in labels):
                        continue
                    f = 1 + P**(parity*n//2)*X**n
                    da = f.gcd(block(v, powers)).degree()
                    assert da > 0
                    if da == n:
                        full += 1
                        continue
                    for z in range(4):
                        y = block(v[:z]+v[z+1:], powers[:z]+powers[z+1:])
                        du = f.gcd(y).degree()
                        weighted = 3*du + 2*da
                        if weighted >= 2*height:
                            high.add((n, v, delta, z, parity, du, da))
                        allocations += 1
    assert (allocations, full) == (64008, 1432)
    assert high == {
        (3, (0, 1, 2, 4), (0, 2, 0, 0), 1, 0, 2, 1),
        (3, (0, 2, 3, 4), (0, 0, 2, 0), 2, 0, 2, 1),
    }
    print(f"{allocations} mixed allocations, {full} full-binomial cases, {len(high)} equality cases")


def large_trinomial_profiles():
    profiles = set()
    for n in range(1, 10):
        for u, v in combinations_with_replacement(range(10), 2):
            if gcd(n, u, v) != 1:
                continue
            for _, s, t in doubled_offsets(P):
                for parity in (0, 1):
                    if parity*n % 2 or (parity*u+s) % 2 or (parity*v+t) % 2:
                        continue
                    y = block((0, u, v), (0, (parity*u+s)//2, (parity*v+t)//2))
                    f = 1+P**(parity*n//2)*X**n
                    mu = f.gcd(y).degree()
                    if 100*mu >= 21*max(n, v):
                        profiles.add((n, u, v, s, t, parity, mu))
    assert profiles == {
        (3, 2, 4, 0, 0, 0, 2), (3, 4, 8, 0, 0, 0, 2),
        (6, 1, 2, 1, 0, 1, 2), (6, 2, 7, 0, 1, 1, 2),
        (6, 5, 7, -1, -1, 1, 2),
        (9, 2, 4, 0, 0, 0, 2), (9, 4, 8, 0, 0, 0, 2),
    }
    maximum = Fraction(0)
    count = 0
    for n, u, v, s, t, parity, mu in sorted(profiles):
        y = block((0, u, v), (0, (parity*u+s)//2, (parity*v+t)//2))
        for r in range(-11, 12):
            h = max(n, v, max(0, r+n)-min(0, r))
            for j in range(-160, 161):
                other = block((0, r, r+n), (0, j, j+parity*n//2))
                dc = y.gcd(other).degree()
                maximum = max(maximum, Fraction(3*mu+2*dc, h))
                count += 1
    assert (count, maximum) == (51681, Fraction(3, 2))
    print(f"{len(profiles)} trinomial profiles, {count} translation checks, maximum ratio {maximum}")


def boundaries():
    h = 10**24
    eps = Fraction(1, 1000)
    q = 2000
    k = 600000
    assert Fraction(1, 48)-Fraction(50, 10**6) > 2*eps
    assert Fraction(2, 3)-2*eps-Fraction(4, h) > 3*Fraction(22, 100)
    assert Fraction(22, 100)-Fraction(6, 10**8)-Fraction(2, h) > Fraction(5, 23)
    assert Fraction(3, 2)+Fraction(18, 10**8)+Fraction(8, h)+2*eps+Fraction(48, 100) < 2
    assert Fraction(5, 3)*(1+Fraction(281, h))+Fraction(18, 10**8)+Fraction(100010, h)+2*eps < 2
    assert 50*(k+1)*q+Fraction(h+2, q) < h*eps
    assert 2*q*(14*(2*k+1)+1) < 10**11
    assert 6+14*(2*k+1)+1+2 < 20000000
    assert Fraction(1, 2*k*32*2000)-Fraction(20000000, h) > Fraction(1, 10**11)
    first_t = 5*10**11
    assert 5*first_t**2 >= h and 20000000 < first_t and 10**11 < 5*first_t
    phi3, phi6 = 1+X+X**2, 1-X+X**2
    assert phi3*phi6 == 1+X**2+X**4 and (1+X)*phi6 == 1+X**3
    assert phi3(-1) == 1 and (X**9-1) % phi3 == 0


if __name__ == "__main__":
    mixed_quadruples()
    large_trinomial_profiles()
    boundaries()
