"""Exclude the endpoint at 1+m*2**k for all k and every five-bit odd m < 2**64."""

from argparse import ArgumentParser
from collections import Counter
from functools import lru_cache
from itertools import combinations, product
from math import comb, gcd, prod

from check_six_binary_digits import endpoint_shapes
from search_even_endpoint import prime_powers


@lru_cache(maxsize=4096)
def divisors(n):
    values = [1]
    for p, w in prime_powers(n).items():
        powers = [1]
        while powers[-1] < w:
            powers.append(powers[-1]*p)
        values = [a*b for a in values for b in powers]
    return values


def normalized_rows(u, factors):
    """All y with sum U and residues (U,0,0) in some order for each full power."""
    t = prod(factors)
    for assignment in product(range(3), repeat=len(factors)):
        if assignment and (assignment[0] or (2 in assignment and
                           (1 not in assignment or assignment.index(2) < assignment.index(1)))):
            continue
        bins = [prod(w for w, j in zip(factors, assignment) if j == i) for i in range(3)]
        least = []
        for group in bins:
            other = t // group
            least.append(other*((u*pow(other, -1, group)) % group) or t)
        room = u - sum(least)
        for x in range(0, room+1, t):
            for y in range(0, room-x+1, t):
                yield (least[0]+x, least[1]+y, least[2]+room-x-y)


def period_rejects(m, a, b, c, u, row, s, t):
    # For this finite check the only shape reaching here has p=23 in T.
    p = 23
    if t % p:
        return False
    modulus = u*s
    value = 2 % modulus
    length = 1
    while value != 1:
        value = 2*value % modulus
        length += 1
    allowed = [k for k in range(length) if (c*pow(2, k, modulus)+1) % modulus == 0]
    # V mod p^2 depends on Q mod p since p divides m.
    for e0 in allowed:
        step = pow(2, length, p)
        q = pow(2, e0, p)
        first = q
        while True:
            v = (m*q+1)*pow(u, -1, p*p) % (p*p)
            if sum((v*y) % (p*p) for y in row) < p*p:
                return False
            q = q*step % p
            if q == first:
                break
    return True


def verify_shapes(m, a, b, c, upper_u, upper_s, factors, counts):
    for u in divisors(upper_u):
        if u < 3 or gcd(u, m) != 1 or 27*(m//upper_s)**2 > u**3:
            continue
        eligible = [p for p, w in factors.items() if a % w == 0 or b % w == 0]
        for mask in range(1 << len(eligible)):
            s = prod(factors[p] for j, p in enumerate(eligible) if mask & (1 << j))
            t = m // s
            t_powers = {p: w for p, w in factors.items() if t % p == 0}
            if 27*t*t > u**3 or any(2*w+1 > u for w in t_powers.values()):
                continue
            for row in normalized_rows(u, tuple(t_powers.values())):
                counts["normalized row pairs"] += 1
                for p, w in prime_powers(u).items():
                    q_digit = -pow(c, -1, p)
                    if sum((x//w * q_digit) % p for x in (a, b)) >= p:
                        counts["carry at a prime of U"] += 1
                        break
                else:
                    for p, w in t_powers.items():
                        v_digit = pow(u, -1, p)
                        if sum((y//w * v_digit) % p for y in row if y % w == 0) >= p:
                            counts["carry at a prime of T"] += 1
                            break
                    else:
                        if not period_rejects(m, a, b, c, u, row, s, t):
                            raise ArithmeticError(("unresolved ray shape", m, a, b, c, u, row, s, t))
                        counts["exponent orbit"] += 1


def check(bits):
    counts = Counter()
    for exponents in combinations(range(1, bits), 4):
        units = (1,) + tuple(1 << i for i in exponents)
        m = sum(units)
        counts["multipliers"] += 1
        shapes = list(endpoint_shapes(units))
        if shapes:
            factors = prime_powers(m)
            for a, b, c, upper_u, upper_s in shapes:
                counts["preliminary binary rows"] += 1
                verify_shapes(m, a, b, c, upper_u, upper_s, factors, counts)
    assert counts["multipliers"] == comb(bits-1, 4)
    print(dict(counts))


if __name__ == "__main__":
    cli = ArgumentParser(description=__doc__)
    cli.add_argument("bits", type=int, nargs="?", default=64)
    check(cli.parse_args().bits)
