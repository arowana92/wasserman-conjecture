"""Check the two-part endpoint bound using exact coefficients and prime powers."""

from itertools import combinations
from math import comb, gcd, prod


def rows(n):
    for x in range(1, n//3+1):
        for y in range(x, (n-x)//2+1):
            yield (x, y, n-x-y), comb(n, x)*comb(n-x, y)


def prime_powers(n):
    p = 2
    while p*p <= n:
        power = 1
        while n % p == 0:
            n //= p
            power *= p
        if power > 1:
            yield p, power
        p += 1
    if n > 1:
        yield n, n


def allocate(atoms, first, second):
    left = right = 1
    partial = False
    for p, power in atoms:
        removed = gcd(power, gcd(first, second))
        retained = power//removed
        partial |= 1 < retained < power
        if gcd(first, power) <= gcd(second, power):
            left *= retained
        else:
            right *= retained
    return left, right, partial


def check_weight(n):
    total = partial = nontrivial = 0
    table = list(rows(n))
    weight, predecessor = tuple(prime_powers(n)), tuple(prime_powers(n-1))
    for row, first in table:
        for other, second in table:
            u, v, part0 = allocate(weight, first, second)
            s, t, part1 = allocate(predecessor, first, second)
            d, j = gcd(n, gcd(first, second)), gcd(n-1, gcd(first, second))
            if u*v*d != n or s*t*j != n-1:
                raise ArithmeticError(("factor assignment", n, row, other))
            if any(c % v for c in other) or prod(c//v for c in other) % (t*t):
                raise ArithmeticError(("complementary parts", n, other, v, t))
            for a, b in combinations(row, 2):
                x, y = a//u, b//u
                if a % u or b % u or x*y % s or x*y*(a+b-1) % (s*s):
                    raise ArithmeticError(("two-part divisibility", n, row, other, a, b))
                if 729*(n-1)**4 > d**6*j**4*(a*b)**3*(a+b-1):
                    raise ArithmeticError(("bound", n, row, other, a, b))
                total += 1
                partial += part0 or part1
                nontrivial += 729*(n-1)**4 > (a*b)**3*(a+b-1)
    return total, partial, nontrivial


def main():
    counts = list(zip(*(check_weight(n) for n in range(3, 46))))
    print("Pair-of-parts checks:", sum(counts[0]))
    print("Checks with a partially retained prime power:", sum(counts[1]))
    print("Checks forcing a gcd greater than one by size:", sum(counts[2]))


if __name__ == "__main__":
    main()
