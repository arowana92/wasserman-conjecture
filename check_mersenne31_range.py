"""Intersect two independent carry sieves on the residual exponent progressions."""

from argparse import ArgumentParser
from bisect import bisect_left, bisect_right
from collections import Counter
from pathlib import Path

from check_mersenne31_exponents import R, factors, free


def prefix_residues(row, p, limit, budget):
    d, u, i, j, a, b, e, h = row
    assert i == 1 and j == 2 and p > 2 and d % p == 0 and pow(2, h, p) == 1
    # N/U is fixed modulo this power of p throughout e + h*t.
    fixed = p
    while pow(2, h, fixed * p) == 1:
        fixed *= p
    mult = d * R
    v = mult
    while v % p == 0:
        fixed *= p
        v //= p
    depth, m = 0, 1
    while m <= limit:
        depth += 1
        m *= p
    precise = fixed * p**depth
    z0 = pow(2, e, precise)
    direction = (R * a, R * b, u - R * (a + b))
    inv = pow(u, -1, precise)
    coeffs = tuple(f * inv % precise for f in direction)
    start_parts = tuple((1 + mult * z0) * f for f in coeffs)
    if not free(start_parts, p, fixed):
        return 1, []
    pairs = [(0, z0)]
    m, w = 1, fixed
    for _ in range(depth):
        step = pow(2, h * m, precise)
        next_w = p * w
        c0, c1, c2 = (f % next_w for f in coeffs)
        updated = []
        for t, z in pairs:
            for digit in range(p):
                n = (1 + mult * z) % next_w
                if (n * c0 % next_w) + (n * c1 % next_w) + (n * c2 % next_w) < next_w:
                    updated.append((t + digit * m, z))
                z = z * step % precise
            if len(updated) > budget:
                break
        if len(updated) > budget:
            break
        pairs = updated
        w = next_w
        m *= p
        if not pairs:
            break
    return m, [t for t, _ in pairs]


def intersection(first, second, limit):
    a, s = first
    b, t = second
    if not s or not t:
        return
    inv = pow(a, -1, b)
    ordered = sorted(v * inv % b for v in t)
    for x in s:
        if x > limit:
            continue
        shift = x * inv % b
        stop = shift + (limit - x) // a
        for k in range(shift // b, stop // b + 1):
            begin = bisect_left(ordered, shift - k * b)
            end = bisect_right(ordered, stop - k * b)
            for y in ordered[begin:end]:
                yield x + a * (y + k * b - shift)


def power_table(e, h, modulus, limit):
    blocks = []
    step = pow(2, h, modulus)
    for _ in range((limit.bit_length() + 7) // 8):
        values = [1]
        for _ in range(255):
            values.append(values[-1] * step % modulus)
        blocks.append(values)
        step = values[-1] * step % modulus
    start = pow(2, e, modulus)

    def evaluate(t):
        result = start
        for values in blocks:
            result = result * values[t & 255] % modulus
            t >>= 8
            if t == 0:
                break
        return result
    return evaluate


def verify(row, bound, budget):
    d, u, i, j, a, b, e, h = row
    if e > bound:
        return 0
    limit = (bound - e) // h
    bases = [p for p, _ in factors(d)]
    assert len(bases) == 2
    if budget is None:
        budget = 850000 if d == 493 else 140000 if d == 481 else 30000
    first, second = [prefix_residues(row, p, limit, budget) for p in bases]
    direction = (a * R, b * R, u - (a + b) * R)
    layers = []
    for digits in (25, 80):
        for p in bases:
            w = p**digits
            inv = pow(u, -1, w)
            layers.append((p, w, inv, power_table(e, h, w, limit)))
    checked = 0
    for t in intersection(first, second, limit):
        checked += 1
        exponent = e + h * t
        for p, modulus, inv, power in layers:
            v = (1 + d * R * power(t)) * inv % modulus
            if not free([z * v for z in direction], p, modulus):
                break
        else:
            raise ArithmeticError(("requires further checks", row, exponent))
    return checked


def cross_checks(row):
    d, u, _, _, a, b, e, h = row
    for p, _ in factors(d):
        modulus, residues = prefix_residues(row, p, 1024, 4096)
        levels = p
        while pow(2, h, p * levels) == 1:
            levels *= p
        temp = d
        while temp % p == 0:
            levels *= p
            temp //= p
        w = levels * modulus
        base = pow(u, -1, w)
        direction = (R * a, R * b, u - R * (a + b))
        allowed = set(residues)
        table = power_table(e, h, w, 1024)
        for t in range(1025):
            z = pow(2, e + h * t, w)
            assert table(t) == z
            scale = (1 + d * R * z) * base % w
            assert free([scale * f for f in direction], p, w) == (t % modulus in allowed)
    first, second = (27, [0, 4, 15, 23]), (25, [1, 2, 12])
    for limit in (0, 50, 300, 4000):
        assert sorted(intersection(first, second, limit)) == [
            t for t in range(limit + 1) if t % 27 in first[1] and t % 25 in second[1]
        ]


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--bound", type=int, default=10**30)
    parser.add_argument("--budget", type=int)
    parser.add_argument("--classes", type=Path, default=Path("mersenne31_exponent_classes.txt"))
    args = parser.parse_args()
    assert args.bound > 0 and (args.budget is None or args.budget > 0)
    totals = Counter()
    for line in args.classes.read_text().splitlines():
        if line.strip() and not line.startswith("#"):
            row = tuple(map(int, line.split()))
            cross_checks(row)
            count = verify(row, args.bound, args.budget)
            totals[row[0]] += count
            print(row[0], row[-2], count, flush=True)
    print("Final candidates excluded by factorial layers:", totals)
