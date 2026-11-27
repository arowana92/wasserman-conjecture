"""Independent checks of primitive-direction recovery and its divisor bounds."""

from functools import lru_cache
from itertools import combinations, combinations_with_replacement, product
from math import comb, gcd, lcm, prod


def profile(x, y):
    n = lcm(sum(x), sum(y))
    H = lcm(*(gcd(*(a*b for a, b in combinations(z, 2))) for z in (x, y)))
    while gcd(n, H) > 1:
        H //= gcd(n, H)
    k = -H % n
    weight = 1 + H // k if gcd(sum(x), sum(y)) == 1 and H % k == 0 else None
    return n, H, weight


@lru_cache(maxsize=None)
def powers(n):
    result = []
    d = 2
    while d*d <= n:
        w = 1
        while n % d == 0:
            w *= d
            n //= d
        if w > 1:
            result.append((d, w))
        d += 1
    if n > 1:
        result.append((n, n))
    return result


def direct_layers(n, rows):
    for offset in (0, 1):
        for _, w in powers(n-offset):
            if all(sum(a % w for a in row) != offset for row in rows):
                return False
    return True


@lru_cache(maxsize=None)
def coefficient(parts):
    n = sum(parts)
    value = 1
    for a in parts:
        value *= comb(n, a)
        n -= a
    return value


def scaled(z, n):
    return tuple(a*(n//sum(z)) for a in z)


def counts():
    rows = [z for length in (2, 3, 4)
            for z in combinations_with_replacement(range(1, 15), length)
            if sum(z) <= 15 and gcd(*z) == 1]
    layers = exact = two = 0
    for x, y in combinations_with_replacement(rows, 2):
        n, H, weight = profile(x, y)
        c = gcd(sum(x), sum(y))
        history = []
        for g in range(1, 6):
            N = g*n
            pair = scaled(x, N), scaled(y, N)
            assert direct_layers(N, pair) == (N == weight)
            layers += 1
            a, b = (coefficient(z) for z in pair)
            G = gcd(N*(N-1), gcd(a, b))
            J = gcd(N-1, gcd(a, b))
            assert G % (c * ((N-1)//gcd(N-1, H))) == 0
            assert H % ((N-1)//J) == 0
            assert N == weight or G > 1
            exact += 1
            for h, old in history:
                assert H * G * old * (g-h) >= c*c*(N-1)*(h*n-1)
                if len(x) == len(y) == 3:
                    assert (n*G*old)**2 > 729*sum(x)*sum(y)*(n*h-1)**2
                two += 1
            history.append((g, G))
    return layers, exact, two


def controls():
    for g, u in product((1, 2, 3, 17, 128, 10**8, 2**127-1, 10**100), (3, 5, 10)):
        h = g*u
        d = h+1
        e = h*d+1
        f = h*d*e-1
        x, y = (1, 1, u-2), (e*f, d*f, d*e)
        n, H, weight = profile(x, y)
        assert gcd(*x) == gcd(*y) == 1 and n == u*sum(y)
        assert weight == g*n and H == weight-1
    x, y = (5, 18), (1, 1)
    assert profile(x, y) == (46, 45, 46)
    assert gcd(46*45, gcd(comb(46, 10), comb(46, 23))) == 1


def units(n, p, limit=3):
    result = []
    place = 1
    while n and len(result) < limit:
        n, digit = divmod(n, p)
        result.extend([place]*min(digit, limit-len(result)))
        place *= p
    return result


def primitive_digit_checks():
    tests = 0
    for n in range(3, 92):
        for x in combinations_with_replacement(range(1, n), 3):
            if sum(x) != n:
                continue
            row_gcd = gcd(*x)
            primitive = tuple(a//row_gcd for a in x)
            value = coefficient(x)
            for p, w in powers(n) + powers(n-1):
                if value % p == 0:
                    continue
                low = units(n, p)
                assert len(low) == 3
                needed = prod(low)
                if n % p == 0:
                    assert row_gcd % w == 0
                    needed //= w**3
                assert prod(primitive) % needed == 0
                if p == 2:
                    assert sum(primitive) % 2 and 32*prod(primitive) < sum(primitive)**3
                tests += 1
    return tests


if __name__ == "__main__":
    controls()
    print("Lowest-layer comparisons, exact endpoint gcds, two-scale bounds:", counts())
    print("Primitive product divisibilities from exact prime-free coefficients:", primitive_digit_checks())
