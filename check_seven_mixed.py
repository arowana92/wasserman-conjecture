"""Exact log-ratio separation and the mixed-prime certificate at seven units."""

from fractions import Fraction
from itertools import permutations
from math import gcd, isqrt


LIMIT = 10**14
PRECISION = 1 << 256


def log_series(x):
    t = (x-1)/(x+1)
    term = t
    value = Fraction(0)
    for k in range(100):
        value += 2*term/(2*k+1)
        term *= t*t
    return value, value + 2*term/(201*(1-t*t))


def log_intervals(primes):
    two_lo, two_hi = log_series(Fraction(2))
    data = {}
    for p in primes:
        shift = p.bit_length()-1
        lo, hi = log_series(Fraction(p, 2**shift))
        lo = (lo+shift*two_lo)*PRECISION
        hi = (hi+shift*two_hi)*PRECISION
        data[p] = lo.numerator//lo.denominator, -(-hi.numerator//hi.denominator)
    return data


def separating_denominator(log_p, log_q):
    low = Fraction(log_p[0], log_q[1])
    high = Fraction(log_p[1], log_q[0])
    mid = (low+high)/2
    n, d = mid.numerator, mid.denominator
    prev_a, a, prev_b, b = 0, 1, 1, 0
    while b <= 2*LIMIT:
        term, rem = divmod(n, d)
        prev_a, a = a, term*a + prev_a
        prev_b, b = b, term*b + prev_b
        n, d = d, rem
    assert gcd(a, b) == 1 and b < 2**70
    assert max(abs(b*low-a), abs(b*high-a))*2*LIMIT < 1
    return b


def valuation_cap(p, q):
    for t in range(1, q):
        if pow(p, t, q) == q-1:
            value = p**t+1
            power = 1
            while value % q == 0:
                power *= q
                value //= q
            return power
    return 0


def solve():
    primes = [p for p in range(3, 256, 2) if all(p%r for r in range(2, isqrt(p)+1))]
    logs = log_intervals(primes)
    solutions = set()
    largest_cap = largest_denominator = 0
    for p, q in permutations(primes, 2):
        cap = valuation_cap(p, q)
        largest_cap = max(largest_cap, cap)
        largest_denominator = max(largest_denominator, separating_denominator(logs[p], logs[q]))
        for u in range(1, 80):
            n = p**u+1
            if n % q:
                continue
            value = n
            j = 0
            while value % q == 0:
                value //= q
                j += 1
            assert cap and q**j <= cap*u
            value -= 1
            e = 0
            while value > 1 and value % q == 0:
                value //= q
                e += 1
            if value == 1 and e:
                solutions.add((p, q, u, j, j+e))
    assert largest_cap == 28561 and 2*largest_cap*80*2**71 < 3**80
    pairs = {(p, q) for p, q, _, _, _ in solutions}
    assert all((q, p) not in pairs for p, q in pairs)
    # Independently construct all sums of two powers q^j + q^k in the same range.
    cross = set()
    by_value = {p**u+1: (p, u) for p in primes for u in range(1, 80)}
    maximum = max(by_value)
    for q in primes:
        small = q
        j = 1
        while small*(q+1) <= maximum:
            k, large = j+1, q*small
            while small+large <= maximum:
                if small+large in by_value:
                    p, u = by_value[small+large]
                    if p != q:
                        cross.add((p, q, u, j, k))
                large *= q
                k += 1
            small *= q
            j += 1
    assert solutions == cross
    print(f'{len(primes)*(len(primes)-1)} log-ratio certificates, '
          f'maximum denominator {largest_denominator}, maximum lifting factor {largest_cap}.')
    print('All (p, q, u, j, k) in 1+p^u = q^j+q^k:', sorted(solutions))


if __name__ == '__main__':
    solve()
