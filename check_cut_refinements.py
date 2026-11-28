"""Check sharp cut constants and the prime-power subdivisions by exact arithmetic."""

from math import comb, gcd

from check_affine_endpoint import factors
from check_direction_pair_count import valuation, divisors
from check_direction_predecessor_loss import constant, vp
from check_part_ratio_divisors import minimum_decrement, sharp_constant


def cut_bound(a, b, d):
    c = sharp_constant(a, b)
    for p, full in factors(b):
        c //= gcd(c, full // gcd(full, d))
    return c


def extend_scale(a, b, p, v):
    size = a + b
    block = p
    while (block - 1) % size:
        block *= p
    shift = p
    while shift <= size*v:
        shift *= p
    return v + shift*((block-1)//size), shift


def check_cut_witnesses(height):
    count = 0
    for a in range(1, height+1):
        for b in range(1, height+1):
            if gcd(a, b) > 1:
                continue
            for p, power in factors(a*b):
                _, v = minimum_decrement(a, b, p, power)
                v, shift = extend_scale(a, b, p, v)
                n = (a+b)*v
                for d in divisors(b):
                    bound = cut_bound(a, b, d)
                    if bound % p:
                        continue
                    s = gcd(d, power) if b % p == 0 else shift
                    row = (a*v, s, b*v-s)
                    assert min(row) > 0 and d % gcd(s, b) == 0
                    assert vp(n-1, p) - valuation(n, row, p) == vp(bound, p)
                    count += 1
    return count


def check_exact(limit):
    count = 0
    for n in range(3, limit+1):
        for unsplit in range(1, n-1):
            v = gcd(n, unsplit)
            a, b = unsplit//v, (n-unsplit)//v
            outer = comb(n, unsplit)
            for s in range(1, (n-unsplit)//2+1):
                coefficient = outer * comb(n-unsplit, s)
                d = gcd(s, b)
                assert coefficient * cut_bound(a, b, d) % (n-1) == 0
                count += 1
    return count


def check_rays():
    transfers = geometric = 0
    for p in (2, 3, 5, 7):
        for m in range(1, 5):
            power = p**m
            for b in range(1, 13):
                for c in range(b, 13):
                    if gcd(b, c) == 1 and b*c*(b+c) % p:
                        assert constant((b+c, power*b, power*c)) == 1
                        transfers += 1
            a = (power-1)//(p-1)
            for x in range(1, power//2+1):
                row = (a, x, power-x)
                upper = gcd(a, x*(power-x))
                if sum(row) <= 350:
                    assert upper % constant(row) == 0
                for v in (1, 2, 3, 5, 12, 31, 101):
                    n = sum(row)*v
                    for q, full in factors(n-1):
                        e = valuation(n, tuple(z*v for z in row), q)
                        assert full <= q**e*gcd(upper, full)
                geometric += 1
    return transfers, geometric


if __name__ == "__main__":
    print("exact coefficient checks", check_exact(85), flush=True)
    print("sharp prime-power witnesses", check_cut_witnesses(20), flush=True)
    print("transfer and geometric directions", *check_rays())
