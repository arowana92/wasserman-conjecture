"""Cofactor 35. Run with plain Python 3.9 or later."""

from functools import lru_cache
from math import gcd

D = 35
PRIMES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31)


def splits():
    allowed = set()
    for i in range(1, D):
        for j in range(i, D):
            for k in range(1, D):
                for l in range(k, (D - 1) // (i + j) - k + 1):
                    if (i * j * k * l) % D == 0 and (i + j) // gcd(i, j) % 2:
                        allowed.add((i, j, k, l))
    assert allowed == {(7, 10, 1, 1)}
    return allowed


def low_classes(bits):
    modulus = 1 << bits
    masks = set()
    for u in range(1, modulus, 2):
        a, b = 7 * u % modulus, 10 * u % modulus
        if a & b:
            continue
        fixed = a | b
        free = modulus - 1 - fixed
        subset = free
        while True:
            masks.add(fixed | subset)
            if subset == 0:
                break
            subset = (subset - 1) & free
    return {q for q in range(1, modulus, 2) if D * q % modulus in masks}


def high_classes(depth):
    lower = 1 << depth
    out = set()
    for w in range(4 * lower, (140 * lower + 16) // 17):
        a, b = w // 10, w // 7
        if a & b:
            continue
        for h in range(max(lower, 17 * w // 70), min(2 * lower, (w + 4) // 4)):
            if (a | b) & ~h == 0:
                out.add(h)
    return out


def density():
    assert low_classes(3) == {5}
    for depth in range(2, 10):
        modulus = 1 << depth
        brute = {q for q in range(1, modulus, 2) if any(
            7 * u % modulus + 10 * u % modulus +
            ((D * q - 17 * u) % modulus) == D * q % modulus and
            not (7 * u % modulus) & (10 * u % modulus) and
            not ((17 * u % modulus) & ((D * q - 17 * u) % modulus))
            for u in range(1, modulus, 2))}
        assert brute == low_classes(depth)
        # Reverse the bounds on w for each h, then enumerate bit positions.
        high = set()
        for h in range(modulus, 2 * modulus):
            for w in range(max(4 * modulus, 4 * h),
                           min((140 * modulus + 16) // 17, (70 * (h + 1) + 16) // 17)):
                a, b = w // 10, w // 7
                if all(((a >> i) & 1) + ((b >> i) & 1) <= ((h >> i) & 1)
                       for i in range(depth + 1)):
                    high.add(h)
                    break
        assert high == high_classes(depth)
    depth = 16
    lower = 1 << depth
    low = sum(q % 4 == 1 for q in low_classes(depth))
    high = high_classes(depth)
    assert (low, len(high)) == (2743, 4243)
    den = 2**30
    upper_sum = sum((den + h - 1) // h for h in high)
    # ln 2 > 69/100 and ln(1+1/h) < 1/h. Bound the product by 0.0085.
    assert low * upper_sum * 100 * 10000 < (lower // 4) * den * 69 * 85
    print("Low classes, high prefixes, upper sum numerator:", low, len(high), upper_sum)


def clear(parts, prime, depth, method):
    if method == "quotients":
        modulus = 1
        for _ in range(depth):
            modulus *= prime
            if sum(v // modulus for v in parts) != sum(parts) // modulus:
                return False
        return True
    work = parts
    for _ in range(depth):
        if sum(v % prime for v in work) >= prime:
            return False
        work = tuple(v // prime for v in work)
    return True


@lru_cache(None)
def scales(p, depth, exponent, offset, allow_nonunit):
    modulus = p**depth
    return sorted({2**offset * pow(z, exponent, modulus) % modulus for z in range(modulus)
                   if z % p or allow_nonunit})


def feasible(x, u, triple, p, depth, exponent, offset, method):
    m = p**depth
    for w in scales(p, depth, exponent, offset, p in (5, 7) and x % p == 0):
        n = 1 + D * w
        cut = (x * w, (D - x) * w + 1)
        if u % p == 0:
            if gcd(n, m) == gcd(u, m) and clear(cut, p, depth, method):
                return True
        else:
            v = n * pow(u, -1, m) % m
            if p in (5, 7) and clear(cut, p, depth, method):
                return True
            if clear(tuple(v * t for t in triple), p, depth, method):
                return True
    return False


def grouped(odd_input, exponent, method):
    # Binary outer power uses Q a square, so the parity of a fixes the
    # square coset of W at every odd prime. For odd t test W=t**exponent.
    rows = set()
    for x in range(3, D):
        if not odd_input and not any(clear((x * z*z, (D-x) * z*z), 2, 4, method)
                                     for z in range(1, 16, 2)):
            continue
        for u in range(3, x + 1):
            if x % u or gcd(u, D) != 1 or (not odd_input and u % 2 == 0):
                continue
            for a in range(1, u // 3 + 1):
                for b in range(a, (u - a) // 2 + 1):
                    f = (a, b, u - a - b)
                    for offset in (range(1) if odd_input else range(2)):
                        if all(feasible(x, u, f, p, 8 if p == 2 else 3, exponent,
                                        offset, method)
                               for p in PRIMES if (p == 2 and odd_input) or (p > 2 and D * u % p == 0)):
                            rows.add((x, u, f, offset))
    return rows


if __name__ == "__main__":
    splits()
    fourths = {pow(z, 4, 27) for z in range(27) if z % 3}
    assert {w for w in fourths if clear((16 * w, 19 * w), 3, 3, "quotients")} == {10, 13, 19, 22}
    density()
    for odd, e in ((False, 2), (True, 4), (True, 8)):
        rows = grouped(odd, e, "quotients")
        assert grouped(odd, e, "digits") == rows
        print("Odd scale, exponent, remaining assignments:", odd, e, len(rows))
        if not odd:
            assert len(rows) == 18
            assert {(x, u, o) for x, u, f, o in rows} == {(9, 3, 0), (9, 9, 0), (11, 11, 0)}
            assert all(not clear((x, D - x), 2, 4, "digits") for x in (9, 11))
        elif e == 4:
            assert rows == {(16, 4, (1, 1, 2), 0)}
        else:
            assert not rows
