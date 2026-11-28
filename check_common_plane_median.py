"""Exhaust the common planes with median absolute coefficient at most 211 above 10^9."""

from collections import Counter
from itertools import combinations_with_replacement
from math import comb, gcd, isqrt, lcm


def prime_factors(n):
    p = 2
    while p * p <= n:
        if n % p == 0:
            yield p
            while n % p == 0:
                n //= p
        p = 3 if p == 2 else p + 2
    if n > 1:
        yield n


def carries(row, p):
    row = list(row)
    while any(row):
        if sum(x % p for x in row) >= p:
            return True
        row = [x // p for x in row]
    return False


def common_prime(n, rows):
    for endpoint in (n, n - 1):
        for p in prime_factors(endpoint):
            if all(carries(row, p) for row in rows):
                return p
    return None


def recover(n, cut, r, s, t):
    first, rem = divmod((t + s) * cut - s * n, r - s)
    second = n - cut - first
    if rem or min(first, second) <= 0:
        return None
    return first, second, cut


def cut_pairs(n, m, t, limit):
    product = m * n
    if limit**2 < product:
        return
    # The sum is 1 + m modulo t. Bound it using the larger root.
    least = isqrt(4 * product - 1) + 1
    greatest = limit + product // limit
    start = least + (1 + m - least) % t
    for total in range(start, greatest + 1, t):
        disc = total**2 - 4 * product
        root = isqrt(disc)
        if root**2 == disc and (total - root) % 2 == 0:
            yield (total - root) // 2, (total + root) // 2


def check_divisors(limit=32):
    count = 0
    for n in range(3, limit + 1):
        rows = []
        for a in range(1, n - 1):
            for x in range(1, (n - a) // 2 + 1):
                y = n - a - x
                u = gcd(x, y)
                coefficient = comb(n, a) * comb(n - a, x)
                rows.append((a, u, coefficient))
        for (a, u, ka), (b, v, kb) in combinations_with_replacement(rows, 2):
            endpoint = n * (n - 1)
            g = gcd(endpoint, ka, kb)
            # The least possible nonzero selected coefficients on separate planes.
            left, right = u // gcd(u, a), v // gcd(v, b)
            assert g * lcm(a, b, left, right) % endpoint == 0
            # The least common positive rational value the two forms can take.
            cross = lcm(u * b, v * a)
            numerator = cross // gcd(cross, a * b)
            assert g * lcm(a, b, numerator) % endpoint == 0
            assert g * cross % endpoint == 0
            count += 1
    return count


def scan(bound=211, cutoff=10**9):
    # The derivation of h=k=1 assumes this inequality.
    assert bound**4 * (cutoff + 1) < 2 * cutoff**2
    counts = Counter()
    for a, b in combinations_with_replacement(range(1, bound + 1), 2):
        if gcd(a, b) != 1 or a == b:
            continue
        m = a * b
        for r, s in ((a, b), (-a, b), (a, -b)):
            positive = max(r, s)
            for t in range((cutoff - 1) // m + 1, positive**2 + 1):
                counts['signed triples in interval'] += 1
                if gcd(t, m) != 1:
                    continue
                n = 1 + m * t
                # Multiply the two integral plane equations after stripping U and V.
                if ((t + s)**2 * m - s**2 * n) % (r - s):
                    continue
                counts['after plane congruence'] += 1
                limit = positive * n // (t + positive)
                for u, v in cut_pairs(n, m, t, limit):
                    counts['integral cut pairs'] += 1
                    if gcd(u, v) != 1:
                        continue
                    rows = recover(n, u, r, s, t), recover(n, v, r, s, t)
                    if None in rows:
                        continue
                    counts['positive row pairs'] += 1
                    if common_prime(n, rows) is None:
                        raise ArithmeticError((n, rows, (r, s, -t)))
    return counts


if __name__ == '__main__':
    print('exact small-weight cut comparisons', check_divisors())
    for key, value in scan().items():
        print(key, value)
