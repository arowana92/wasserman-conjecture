"""Direct checks of the forced row near a leading predecessor power."""

from math import comb, gcd


def prime_power_divisors(n):
    factor = 2
    while factor * factor <= n:
        power = 1
        while n % factor == 0:
            n //= factor
            power *= factor
            yield factor, power
        factor += 1
    if n > 1:
        yield n, n


def adds_freely(parts, prime):
    total = sum(parts)
    while total:
        if sum(v % prime for v in parts) != total % prime:
            return False
        total //= prime
        parts = tuple(v // prime for v in parts)
    return True


def coefficient(parts):
    a, b, c = parts
    return comb(a + b + c, a) * comb(b + c, b)


def main():
    assignments = matching = clear = 0
    for p in (2, 3, 5, 7, 11):
        for m in range(3, 241):
            high = 1
            while p * high <= m:
                high *= p
            powers = [(q, w) for q, w in prime_power_divisors(m)
                      if q != p and m - high < 2 * w]
            if not powers:
                continue
            for x in range(1, m // 2 + 1):
                for y in range(x, m - x):
                    z = m - x - y
                    g = gcd(x, y)
                    if g < 3 or not adds_freely((x, y, z), p):
                        continue
                    for a in range(1, 6):
                        shift = p**a
                        u = gcd(g, 1 + shift * m)
                        if u < 3:
                            continue
                        parts = (shift * x, shift * y, shift * z + 1)
                        for q, w in powers:
                            assignments += 1
                            assert u < 2 * w
                            if sorted(v % w for v in parts) != [0, 0, 1]:
                                continue
                            matching += 1
                            r = pow(shift, -1, w)
                            forced = (shift * (m - w), shift * r, shift * (w - r) + 1)
                            assert sorted(parts) == sorted(forced)
                            assert gcd(m // w - 1, r) % u == 0
                            if adds_freely(parts, q):
                                clear += 1
                                assert gcd(coefficient(parts), p * q) == 1
    assert (assignments, matching, clear) == (47774, 12, 4)
    assert gcd(573 * 572, coefficient((12, 528, 33)), coefficient((191,) * 3)) == 39
    print("Assignments, matching rows, rows avoiding both primes:", assignments, matching, clear)


if __name__ == "__main__":
    main()
