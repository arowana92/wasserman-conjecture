"""Enumerate the exponent classes left at d*(2^31-1)*2^e+1 for d<=512."""

from argparse import ArgumentParser
from collections import Counter
from functools import lru_cache
from itertools import combinations
from math import gcd, isqrt, lcm
from pathlib import Path
from subprocess import run

from check_binary_gap_finite import factors
from check_mersenne_multipliers import shapes

R = 2**31 - 1


@lru_cache(None)
def order(modulus):
    phi = modulus
    for p, _ in factors(modulus):
        phi = phi // p * (p - 1)
    h = phi
    for p, _ in factors(phi):
        while h % p == 0 and pow(2, h // p, modulus) == 1:
            h //= p
    return h


def merge(a, h, b, k):
    g = gcd(h, k)
    if (b - a) % g:
        return None
    n = k // g
    t = (b - a) // g * pow(h // g, -1, n) % n
    return a + h * t, h * n


def subgroup_log(base, target, size, modulus):
    step = isqrt(size) + 1
    digits = {}
    value = 1
    for j in range(step):
        if value == target:
            return j
        digits[value] = j
        value = value * base % modulus
    backward = pow(value, -1, modulus)
    value = target
    for i in range(1, step + 1):
        value = value * backward % modulus
        if value in digits and i * step + digits[value] < size:
            return i * step + digits[value]
    return None


@lru_cache(None)
def weight_class(d, u):
    h = order(u)
    target = -pow(d * R, -1, u) % u
    if pow(target, h, u) != 1:
        return None
    e, period = 0, 1
    for _, w in factors(h):
        b = subgroup_log(pow(2, h // w, u), pow(target, h // w, u), w, u)
        if b is None:
            return None
        e, period = merge(e, period, b, w)
    assert period == h
    return (e, h) if pow(2, e, u) == target else None


def free(parts, p, w):
    while w > 1:
        if sum(x % w for x in parts) >= w:
            return False
        w //= p
    return True


def initial_layers(d, u, i, j):
    for p, w in factors(u):
        value = -pow(d * R, -1, w) * (u // w)
        if not free((i * value, j * value), p, w):
            return False
    return True


def passes(shape, p, w, z):
    d, u, i, j, a, b = shape
    t = d * R
    n = t * z + 1
    first = (u * i * z, u * j * z, (t - (i + j) * u) * z + 1)
    if u % p == 0:
        return free(first, p, w)
    if n % p not in (0, 1):
        return True
    scale = n * pow(u, -1, w) % w
    second = (scale * R * a, scale * R * b, scale * (u - R * (a + b)))
    return free(second, p, w) or (n % p == 1 and free(first, p, w))


def local_orbit(shape, e, h, p, cap):
    w = p
    while w * p <= cap:
        w *= p
    u = shape[1]
    if u % p == 0:
        power = next(q for r, q in factors(u) if r == p)
        if w <= power:
            w = power * p
    length = order(w) // gcd(order(w), h)
    if length > 128:
        return True
    value, step = pow(2, e, w), pow(2, h, w)
    for _ in range(length):
        if passes(shape, p, w, value):
            return True
        value = value * step % w
    return False


def deeper_layers(shape, e, h, p):
    u = shape[1]
    period = lcm(h, order(p))
    w = p
    while pow(2, period, w * p) == 1:
        w *= p
    full = next((q for r, q in factors(u) if r == p), 1)

    def viable(z, modulus):
        if full > 1 and modulus > full and (shape[0] * R * z + 1) % (p * full) == 0:
            return False
        return passes(shape, p, modulus, z)

    states = {pow(2, e + t * h, w) for t in range(period // h)}
    states = {z for z in states if viable(z, w)}
    for _ in range(35):
        if not states or len(states) > 150:
            break
        new = w * p
        states = {z + t * w for z in states for t in range(p) if viable(z + t * w, new)}
        w = new
    return bool(states)


def all_classes(shape, mask):
    d, u, i, j, _, _ = shape
    if not initial_layers(d, u, i, j):
        return
    data = weight_class(d, u)
    if data is None:
        return
    for b in range(order(d)):
        z = pow(2, b, d)
        if not all(passes(shape, p, w, z) for p, w in factors(d)):
            continue
        combined = merge(*data, b, order(d))
        if combined is None:
            continue
        for j0 in range(31):
            if mask >> j0 & 1:
                result = merge(*combined, j0, 31)
                if result is not None:
                    yield result


def check_small_scales(reported):
    # The all-height theorem leaves only these six shapes through odd d=127.
    # For i=1,j=2 the binary condition says that U has no adjacent ones.
    cases = {105: [(5, 21), (7, 15), (10, 21), (14, 15)], 117: [(9, 13), (13, 18)]}
    direct = set()
    for d, companions in cases.items():
        t = d * R
        places = t & (t >> 1)
        candidates = [0]
        for bit in range(places.bit_length()):
            if places >> bit & 1:
                power = 1 << bit
                candidates += [u + power for u in candidates if not (2 * u & power)]
        for u in candidates:
            if u % 2 == 0 or 3 * u >= t or gcd(u, d * R) != 1:
                continue
            inverse = pow(u, -1, R)
            for a, b in companions:
                if u <= R * (a + b):
                    continue
                # No odd prime of d divides either normalized first slope.
                if any(sorted((R * a % w, R * b % w, (u - R * (a + b)) % w))
                       != [0, 0, u % w] for _, w in factors(d)):
                    continue
                z = (inverse * a) % R + (inverse * b) % R
                mask = sum(1 << e for e in range(31) if z <= (d << e) % R)
                if mask:
                    direct.add((d, u, 1, 2, a, b, mask))
    assert direct == reported


def rows(checker):
    lines = []
    for d in range(3, 512, 2):
        powers = [w for _, w in factors(d)]
        lines.extend(" ".join(map(str, (d, *shape))) for shape in shapes(d, powers))
    result = run([str(checker.resolve())], input="\n".join(lines) + "\n",
                 text=True, capture_output=True, check=True)
    small = [p for p in range(3, 98, 2) if len(factors(p)) == 1 and factors(p)[0][0] == p]
    counts = Counter()
    residual = []
    small_scales = set()
    for line in result.stdout.splitlines():
        *shape, mask = map(int, line.split())
        d, u = shape[:2]
        if d <= 127:
            small_scales.add((*shape, mask))
        counts["scale and shape choices"] += 1
        forced = [p for p, _ in factors(d)] + [p for p, _ in factors(u) if p <= 257]
        for e, h in all_classes(shape, mask):
            counts["joint exponent classes"] += 1
            if not all(local_orbit(shape, e, h, p, cap)
                       for p in dict.fromkeys(forced + small) for cap in (64,256,1024,4096,16384)):
                continue
            counts["after finite orbits"] += 1
            if all(deeper_layers(shape, e, h, p) for p in forced
                   if p <= 101 or d % p == 0):
                residual.append((*shape, e, h))
    check_small_scales(small_scales)
    return residual, counts


def residue_obstruction(e, h, exponent):
    # An explicit prime-power divisor of h suffices to exclude every kth power.
    for _, w in factors(h):
        for n in range(2, min(w, 4096) + 1):
            if w % n == 0 and e % n not in {pow(v, exponent, n) for v in range(n)}:
                return True
    return False


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--scale-checker", type=Path, required=True)
    parser.add_argument("--compare", type=Path)
    args = parser.parse_args()
    residual, counts = rows(args.scale_checker)
    print(counts)
    for row in residual:
        print(*row)
    for first, second in combinations(residual, 2):
        assert merge(first[-2], first[-1], second[-2], second[-1]) is None
    for first in residual:
        for second in residual:
            assert merge(first[-2] + 1, first[-1], second[-2], second[-1]) is None
    assert all(residue_obstruction(e, h, k) for *_, e, h in residual for k in (2, 3))
    if args.compare:
        expected = [tuple(map(int, line.split())) for line in args.compare.read_text().splitlines()
                    if line.strip() and not line.startswith("#")]
        assert sorted(expected) == sorted(residual)
