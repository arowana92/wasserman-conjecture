"""Check the examples and the digit-sum obstruction for three predecessor primes."""

from collections import deque
from functools import lru_cache
from itertools import combinations
from math import gcd, isqrt, prod


def digit_sum(value, base):
    total = 0
    while value:
        value, digit = divmod(value, base)
        total += digit
    return total


@lru_cache(None)
def minimum(d, p):
    """Shortest nonempty sum of powers of p congruent to zero modulo d."""
    assert d > 1 and gcd(d, p) == 1
    residues = {1}
    x = p % d
    while x != 1:
        residues.add(x)
        x = x * p % d
    # Divide all summands by the first power, making its residue 1.
    queue = deque([(1, 1)])
    seen = {1}
    while queue:
        value, cost = queue.popleft()
        for r in residues:
            nxt = (value + r) % d
            if nxt == 0:
                return cost + 1
            if nxt not in seen:
                seen.add(nxt)
                queue.append((nxt, cost + 1))
    raise ArithmeticError(d)


def prime(n):
    return n > 1 and all(n % i for i in range(2, isqrt(n) + 1))


def examples():
    mersenne = 2**31 - 1
    shifted = mersenne - (2**14 + 2**2) + 2**45 + 2**64
    assert prime(mersenne) and prime(8589950981)
    assert shifted == mersenne * 8589950981 == 18446779260229107707
    assert digit_sum(shifted, 2) == 31
    t = 131071 * 131149
    assert prime(131071) and prime(131149) and 131071 == 2**17 - 1
    assert t == 17189830579
    assert digit_sum(t, 2) == 19 and digit_sum(1755 * t, 2) == 17


def rows_avoid(total, row, p):
    power = p
    while power <= total:
        if total // power != sum(a // power for a in row):
            return False
        power *= p
    return True


def small_checks():
    for p in (2, 3, 5, 7):
        for d in range(2, 46):
            if gcd(p, d) != 1:
                continue
            mu = minimum(d, p)
            # Enumerate base-p digit strings by successive residue transitions.
            # A minimum-cost positive multiple has a representative with at
            # most d digits: delete cycles in the path of nonzero prefixes.
            states = {k % d: k for k in range(1, min(d, p))}
            if p > d:
                states[0] = d
            for _ in range(d - 1):
                updated = dict(states)
                for v, cost in states.items():
                    for digit in range(p):
                        dest = (p * v + digit) % d
                        updated[dest] = min(updated.get(dest, 10**9), cost + digit)
                states = updated
            assert states[0] == mu

    tested = feasible = 0
    for p in (2, 3, 5, 7):
        atoms = [
            (q**j, q)
            for q in (2, 3, 5, 7, 11, 13, 17, 19)
            if q != p
            for j in range(1, 6)
            if q**j <= 50
        ]
        for (a, pa), (b, pb) in combinations(atoms, 2):
            if pa == pb:
                continue
            for exponent in range(1, 8):
                powers = (p**exponent, a, b)
                m = prod(powers)
                row = tuple((m // w) * pow(m // w, -1, w) for w in powers)
                tested += 1
                if sum(row) == m + 1 and rows_avoid(m + 1, row, p):
                    feasible += 1
                    assert digit_sum(a * b, p) + 1 >= minimum(a, p) + minimum(b, p) + minimum(
                        a * b, p
                    )
    assert feasible
    print(tested, "CRT assignments checked,", feasible, "with correct sum and no base-p carry")


if __name__ == "__main__":
    examples()
    small_checks()
    print("Both all-exponent examples verified")
