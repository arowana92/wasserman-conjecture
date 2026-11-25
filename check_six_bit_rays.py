"""Exclude the endpoint at 1+m*2**k for all k and every five-bit odd m < 2**96."""

from argparse import ArgumentParser
from collections import Counter
from functools import lru_cache
from itertools import combinations, product
from math import comb, gcd, prod

from flint import fmpz

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


# For each remaining shape, a prime, number of low digits, and orbit period.
# The verifier tests both the period and the carry at every point of the orbit.
ORBIT_CHECKS = {
    6917529027643179011: (23, 2, 11),
    885443715538058477603: (11, 3, 55),
    3626777458843904703987715: (5, 3, 1),
    3022316854879582150459397: (29, 2, 14),
    15010022976335235833151881219: (139, 2, 3),
}


def period_rejects(m, a, b, c, u, row, s, t):
    if m not in ORBIT_CHECKS:
        return False
    p, digits, period = ORBIT_CHECKS[m]
    w = p**digits
    if (u*s*t) % p or not fmpz(p).is_prime():
        return False
    modulus = u*s
    value = 2 % modulus
    length = 1
    while value != 1:
        value = 2*value % modulus
        length += 1
    allowed = [k for k in range(length) if (c*pow(2, k, modulus)+1) % modulus == 0]
    # For the C row, only mQ mod w is needed to compute V mod w.
    use_c = t % p == 0
    q_modulus = w // gcd(w, m) if use_c else w
    if pow(2, length*period, q_modulus) != 1 % q_modulus:
        return False
    for e0 in allowed:
        for j in range(period):
            q = pow(2, e0+length*j, q_modulus)
            if use_c:
                v = (m*q+1)*pow(u, -1, w) % w
                parts = [v*y for y in row]
            else:
                parts = (q*a, q*b, q*c+1)
            # A completed full power of U also requires V to be a p-unit.
            common = gcd(u, w)
            if 1 < common < w and (m*q+1) % (p*common) == 0:
                continue
            if all(sum(x % (p**i) for x in parts) < p**i for i in range(1, digits+1)):
                return False
    return True


@lru_cache(maxsize=4096)
def prime_digits(u):
    return prime_powers(u)


def verify_shapes(m, a, b, c, upper_u, upper_s, counts, count_raw_pairs):
    for u in divisors(upper_u):
        if u < 3 or gcd(u, m) != 1 or 27*(m//upper_s)**2 > u**3:
            continue
        bad_u = any(sum((-(x//w)*pow(c, -1, p)) % p for x in (a, b)) >= p
                    for p, w in prime_digits(u).items())
        if bad_u and not count_raw_pairs:
            continue
        factors = prime_digits(m)
        eligible = [p for p, w in factors.items() if a % w == 0 or b % w == 0]
        for mask in range(1 << len(eligible)):
            s = prod(factors[p] for j, p in enumerate(eligible) if mask & (1 << j))
            t = m // s
            t_powers = {p: w for p, w in factors.items() if t % p == 0}
            if 27*t*t > u**3 or any(2*w+1 > u for w in t_powers.values()):
                continue
            for row in normalized_rows(u, tuple(t_powers.values())):
                counts["normalized row pairs"] += 1
                if bad_u:
                    counts["carry at a prime of U"] += 1
                    continue
                for p, w in t_powers.items():
                    v_digit = pow(u, -1, p)
                    if sum((y//w * v_digit) % p for y in row if y % w == 0) >= p:
                        counts["carry at a prime of T"] += 1
                        break
                else:
                    if not period_rejects(m, a, b, c, u, row, s, t):
                        raise ArithmeticError(("unresolved ray shape", m, a, b, c, u, row, s, t))
                    counts["exponent orbit"] += 1


def check(bits, count_raw_pairs=False):
    counts = Counter()
    for exponents in combinations(range(1, bits), 4):
        units = (1,) + tuple(1 << i for i in exponents)
        m = sum(units)
        counts["multipliers"] += 1
        for a, b, c, upper_u, upper_s in endpoint_shapes(units):
            counts["preliminary binary rows"] += 1
            verify_shapes(m, a, b, c, upper_u, upper_s, counts, count_raw_pairs)
    assert counts["multipliers"] == comb(bits-1, 4)
    print(dict(counts))


if __name__ == "__main__":
    cli = ArgumentParser(description=__doc__)
    cli.add_argument("bits", type=int, nargs="?", default=96)
    cli.add_argument("--count-raw-pairs", action="store_true",
                     help="also enumerate CRT rows already excluded by a prime of U")
    args = cli.parse_args()
    check(args.bits, args.count_raw_pairs)
