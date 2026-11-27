"""Exact controls for growing_prime_power_steps.md (standard library)."""

from itertools import combinations
from math import gcd, prod

from check_endpoint_pair_bound import prime_powers, rows


def digit_margin(p):
    # Clear denominators in (sum(x**j) * x**(-(p-1)/6))**18 < p**16.
    w, v = p-1, p+5
    return (v**p-w**p)**18 < (v-w)**18 * v**(15*w) * w**(3*w) * p**16


def check_digit_margins():
    primes = [p for p in range(2, 1250) if list(prime_powers(p)) == [(p, p)]]
    assert all(digit_margin(p) for p in primes)
    assert not digit_margin(1259)
    print(f"Exact digit inequality passed at all {len(primes)} primes through 1249.")


def factorial_valuation(n, p):
    value = 0
    while n:
        n //= p
        value += n
    return value


def check_shifted_rows():
    count = shared = 0
    for n in range(3, 43):
        for slopes, coefficient in rows(n):
            for p, scale in ((2, 2), (2, 16), (3, 9), (5, 125), (7, 7)):
                for slot in range(3):
                    parts = [a*scale + (i == slot) for i, a in enumerate(slopes)]
                    valuation = factorial_valuation(sum(parts), p) - sum(
                        factorial_valuation(z, p) for z in parts)
                    small = coefficient
                    expected = 0
                    while small % p == 0:
                        expected += 1
                        small //= p
                    assert valuation == expected
                    count += 1
                    shared += n % p == 0
    print(f"{count:,} shifted factorial valuations, including {shared:,} with p dividing n.")


def reconstruct(scale, s, bins, v, parts):
    total = sum(parts)
    u = total//v
    t_all = prod(bins)
    n = s*t_all
    assert total == scale*n+1 and gcd(u, scale*s) == 1
    # Calculate the residue using only U, s and the fixed scale.
    assert n % (u*s) == s * (-pow(scale*s, -1, u) % u)
    index = min(range(3), key=bins.__getitem__)
    d = bins[index]
    t = t_all//d
    norm = parts[index]//v
    z, rz = divmod(norm, t)
    y, ry = divmod(u-norm, d)
    assert rz == ry == 0 and min(z, y) >= 1
    dividend = scale*s*d*d*y-z
    assert dividend > 0 and dividend == (scale*s*d-v*z)*u
    assert s*d*((u-d*y)//z) == n and (u-d*y) % z == 0


def check_complements():
    exact = built = 0
    for total in range(4, 96):
        pred = tuple(prime_powers(total-1))
        weight = tuple(prime_powers(total))
        for parts, coefficient in rows(total):
            v = prod(w for p, w in weight if coefficient % p)
            if min(v, total//v) < 3:
                continue
            for base, full in pred:
                bins = [1, 1, 1]
                for p, w in pred:
                    if p != base and coefficient % p:
                        residues = [z % w for z in parts]
                        assert sorted(residues) == [0, 0, 1]
                        bins[residues.index(1)] *= w
                scale = base
                while scale <= full:
                    s = (total-1)//(scale*prod(bins))
                    reconstruct(scale, s, bins, v, parts)
                    scale *= base
                    exact += 1
    for bins in combinations(range(2, 16), 3):
        if any(gcd(i, j) > 1 for i, j in combinations(bins, 2)):
            continue
        t = prod(bins)
        for v in range(3, 19):
            if gcd(t, v) > 1:
                continue
            parts = [v*(t//d)*pow(v*(t//d), -1, d) for d in bins]
            for base, scale in prime_powers((sum(parts)-1)//t):
                if t % base:
                    s = (sum(parts)-1)//(t*scale)
                    reconstruct(scale, s, bins, v, parts)
                    built += 1
    assert min(exact, built) > 100
    print(f"{exact:,} complementary rows/scales from exact coefficients and "
          f"{built:,} reconstructions with three nonempty bins.")


if __name__ == "__main__":
    check_digit_margins()
    check_shifted_rows()
    check_complements()
