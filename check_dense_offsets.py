"""Check the dense offset theorems by direct allocation of ordinary base-p digits."""

from collections import Counter
from functools import lru_cache
from itertools import combinations_with_replacement
from math import gcd, prod


@lru_cache(maxsize=None)
def splits(digit):
    return [(a, b, digit-a-b) for a in range(digit+1) for b in range(digit-a+1)]


def free_rows(n, prime):
    rows = [(0, 0, 0)]
    place = 1
    while n:
        n, digit = divmod(n, prime)
        if digit:
            rows = [(a+place*x, b+place*y, c+place*z)
                    for x, y, z in splits(digit) for a, b, c in rows]
        place *= prime
    return (row for row in rows if 0 < row[0] <= row[1] <= row[2])


def has_size_obstruction(n, row):
    # Use all row gcds directly. No choice of the assigned prime factors is needed.
    g = gcd(*row)
    if g < 3:
        return True
    a, b, c = row
    pred = n-1
    bins = (gcd(pred, a, b), gcd(pred, a, c), gcd(pred, b, c))
    return g**3 * prod(bins)**2 < 27*pred**2


def digit_controls():
    from math import factorial

    for n in range(3, 40):
        for p in (2, 3, 5, 7):
            direct = {row for row in combinations_with_replacement(range(1, n-1), 3)
                      if sum(row) == n and factorial(n)//prod(map(factorial, row)) % p}
            assert set(free_rows(n, p)) == direct
    assert not has_size_obstruction(3000, (1000, 1000, 1000))


def sample_weights():
    # Overlapping parameter choices can describe the same (weight, prime).
    weights = set()
    for p in (2, 3, 5, 7):
        for m in (2, 5, 10, 20):
            for j in (2, 3, 4, 5, 6, 7):
                q = p**(j*m)
                for i in range(1, j):
                    middle = p**(i*m)
                    for t in (1, 2, 3, 4, 7, 10, 15):
                        for d in range(1, 33):
                            if d*(d-1) % p:
                                continue
                            allowed = d**(5*j) <= q*q if t == 3 else (
                                j%2 == 0 and i%2 == 1 and d >= 2 and
                                t**(5*j+2)*d**(5*j) <= q*q)
                            if allowed:
                                weights.add((t*q+middle+d, p))
        for k in (51, 77, 103):
            q = p**k
            for lower_exponent in (1, 2, 9, 21, 47, k-2, k-1):
                middle = p**lower_exponent
                for d in range(1, 33):
                    if d**25 <= q and not d*(d-1)%p and d*middle < q and (d<middle or p==2):
                        weights.add((3*q+middle+d, p))
    return weights


def check():
    digit_controls()
    by_prime = Counter()
    checked = 0
    for n, p in sorted(sample_weights()):
        for row in free_rows(n, p):
            assert has_size_obstruction(n, row), (n, p, row)
            checked += 1
        by_prime[p] += 1
    print('Weights by assigned prime:', sorted(by_prime.items()))
    print('Weight and prime pairs:', sum(by_prime.values()))
    print('Positive unordered carry-free rows:', checked)


if __name__ == '__main__':
    check()
