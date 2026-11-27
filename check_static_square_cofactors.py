"""Compare quadratic residue and exact-gcd tests at one or two times a prime power."""

from functools import lru_cache
from math import gcd

from check_static_power_transfer import prime_powers
from check_three_factors import powers


@lru_cache(None)
def factors(x):
    return tuple(powers(x))


def possible_square(h, b):
    # Enumerating the unit squares also treats powers of two.
    squares = {t*t % b for t in range(b) if gcd(t, b) == 1}
    for x in range(3, h*b):
        if x % b not in squares:
            continue
        if x % 4 == 0 and (h*b) % 4 == 3:
            return x, 4
        for ell in factors(x):
            if ell != 2 and gcd(ell, h*b) == 1 and pow(-h*b, (ell-1)//2, ell) == 1:
                return x, ell
    return None


def exact_table(h, b):
    for x in range(3, h*b):
        if any(x*z*z % b == 1 for z in range(b)):
            if any(gcd(x, h*b*t*t+1) > 2 for t in range(x)):
                return False
    return True


def check():
    expected = {1: [2, 3, 4, 5, 8, 9, 13, 16, 25, 37], 2: [3, 5, 9, 11, 29]}
    for h in (1, 2):
        certified = []
        for b, q in sorted(prime_powers(1000).items()):
            if h == 2 and q == 2:
                continue
            excludes = possible_square(h, b) is None
            if b <= 100:
                assert excludes == exact_table(h, b)
            if excludes:
                certified.append(b)
        assert certified == expected[h]
        print(f"square exclusions at multiplier {h}, prime powers through 1000: {certified}")
    assert possible_square(2, 7) == (9, 3)
    for w in range(9*7*5):
        for x in range(3, 10):
            if x*w % 5 != 1:
                continue
            for u in (3, 7, 9):
                n = 10*w+1
                if x % u or n % u or gcd(n//u, u) != 1:
                    continue
                if u == 9:
                    assert x == 9 and (w+1) % 9 == 0
                elif u == 3:
                    assert w % 9 in (2, 5)
                else:
                    assert x == 7 and w % 7 == 2 and w % 5 == 3
    assert {t**3 % 9 for t in range(9)} == {0, 1, 8}
    assert {t**3 % 7 for t in range(7)} == {0, 1, 6}
    print("scale alternatives at ten and cubic residue exclusions verified")


if __name__ == "__main__":
    check()
