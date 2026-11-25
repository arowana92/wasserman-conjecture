"""Certify the finite part of the theorem for N = 5 * 2**k + 1 (Python 3.9+)."""

from bisect import bisect_right
from math import gcd, isqrt, lcm, prod
from pathlib import Path


def certified_primes():
    """Prove the small primes by trial division and the larger ones by order."""
    known = {n for n in range(2, 65536)
             if all(n % d for d in range(2, isqrt(n) + 1))}
    path = Path(__file__).with_name("binary_ray_primes.txt")
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        p, witness, *factors = map(int, line.split())
        remainder = p - 1
        if p < 65536 or p in known or pow(witness, remainder, p) != 1:
            raise ArithmeticError(("invalid order certificate", p))
        for q in factors:
            if q not in known or remainder % q or pow(witness, (p - 1) // q, p) == 1:
                raise ArithmeticError(("invalid order factor", p, q))
            while remainder % q == 0:
                remainder //= q
        if remainder != 1:
            raise ArithmeticError(("incomplete order factorization", p))
        known.add(p)
    return sorted(known)


def offset_constants():
    height = 15
    small_layers = lcm(*range(1, height + 1))
    return {
        h: lcm(small_layers, *(abs(h - 1 - 5*u) for u in range(h + 1)))
        for h in range(height) if h % 5 != 1
    }


def finite_exponents():
    constants = offset_constants()
    degree = sum(15 - h for h in constants)
    if degree != 93:
        raise ArithmeticError(degree)
    bound = (2**degree) * prod(kh**(15-h) for h, kh in constants.items())
    # For any larger k, N**3 > bound. This also holds for some smaller k.
    cutoff = (bound.bit_length() + 2) // 3
    remaining = []
    for k in range(1, cutoff):
        n = (5 << k) + 1
        if n < 30:
            remaining.append(k)
            continue
        divisor = prod(((n - h) // gcd(n - h, kh))**(15-h)
                       for h, kh in constants.items())
        if divisor <= n**90:
            remaining.append(k)
    return cutoff, remaining


def has_carry(row, p):
    a, b, c = row
    while a or b or c:
        if a % p + b % p + c % p >= p:
            return True
        a, b, c = a // p, b // p, c // p
    return False


def factorial_valuation(n, p):
    result = 0
    while n:
        n //= p
        result += n
    return result


def valuation(n, row, p):
    return factorial_valuation(n, p) - sum(factorial_valuation(x, p) for x in row)


def check_remaining(k, primes):
    t = 2**k
    n = 5*t + 1
    r = primes[bisect_right(primes, n) - 1]
    if r <= 4*t:
        raise ArithmeticError(("missing nearby prime", k))
    gap = n - r
    binary_row = (1, t, 4*t)
    checked = 0
    for i in range(1, gap // 2 + 1):
        for j in range(i, gap - i + 1):
            near_row = (i, j, n - i - j)
            u = gcd(n, i, j)
            v = n // u
            for a in range(1, u // 3 + 1):
                for b in range(a, (u - a) // 2 + 1):
                    last_row = (v*a, v*b, v*(u-a-b))
                    rows = (binary_row, near_row, last_row)
                    common = {p for p in (3, 5) if all(has_carry(row, p) for row in rows)}
                    by_valuation = {p for p in (3, 5) if all(valuation(n, row, p) > 0 for row in rows)}
                    if not common or common != by_valuation:
                        raise ArithmeticError((n, rows, common, by_valuation))
                    checked += 1
    return checked


def main():
    primes = certified_primes()
    cutoff, remaining = finite_exponents()
    tested = sum(check_remaining(k, primes) for k in remaining)
    print(f"The uniform bound excludes k >= {cutoff}.")
    print(f"Exact products leave {len(remaining)} positive exponents, largest {max(remaining)}.")
    print(f"Verified {tested} candidate row pairs by carries and factorial valuations.")
    print("For k = 0, N = 6 has only two binary digit units and 2 is universal.")


if __name__ == "__main__":
    main()
