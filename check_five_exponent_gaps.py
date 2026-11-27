"""Prove exclusions for gaps between binary exponents by quinary carries."""
import argparse
import json
from pathlib import Path

BANDS = ((0, 6), (10, 12), (15, 18), (20, 24))  # denominator 30
GUARD = 7
LIMIT = 4000


def low_prefixes(multiplier, depth, cap=None):
    modulus = 5
    front = [2]
    for column in range(2, depth + 1):
        old = modulus
        modulus *= 5
        ratio = multiplier % modulus
        following = []
        for prefix in front:
            for x in range(prefix, modulus, old):
                if 3*x % modulus + (2*x + 1) % modulus >= modulus:
                    continue
                y = x*ratio % modulus
                if 3*y % modulus + (2*y + 1) % modulus < modulus:
                    following.append(x)
        front = following
        if not front or (cap is not None and len(front) > cap):
            return column, front
    return depth, front


def carry_free(t):
    a, b = 3*t, 2*t + 1
    while a or b:
        if a % 5 + b % 5 >= 5:
            return False
        a //= 5
        b //= 5
    return True


def clip_bands(lo, hi, step, pad):
    for q in range((lo - pad) // (30*step), hi // (30*step) + 1):
        for a, b in BANDS:
            left = max(lo, (30*q + a)*step)
            right = min(hi, (30*q + b)*step + pad)
            if left <= right:
                yield left, right


def high_prefixes(multiplier, height, cap=2200):
    scale = 5**height
    near = 1
    while 5*near <= multiplier:
        near *= 5
    # Integer coordinates are 30*multiplier*scale times the first mantissa.
    front = [(6*multiplier*scale, 30*multiplier*scale)]
    pad = (15*multiplier + 5**GUARD - 1) // 5**GUARD
    step = 5*scale
    for level in range(-1, height + 1):
        following = []
        for lo, hi in front:
            for left, right in clip_bands(lo, hi, multiplier*step, 0):
                following.extend(clip_bands(left, right, near*step, pad))
        front = following
        if not front or len(front) > cap:
            return level, len(front)
        step //= 5
    return height, len(front)


def gap_certificates():
    proofs = {}
    for gap in range(4, LIMIT + 1, 4):
        multiplier = 1 << gap
        depth, residues = low_prefixes(multiplier, 262 if gap == 36 else 105,
                                       2000 if gap == 36 else 300)
        if not residues:
            proofs[gap] = ('low', depth)
            continue
        height, size = high_prefixes(multiplier, 320)
        if size:
            proofs[gap] = ('open',)
            continue
        # Recompute with the shorter scale and its own perturbation bound.
        assert not high_prefixes(multiplier, height)[1]
        _, short = low_prefixes(multiplier, height + GUARD)
        assert all(not (carry_free(t) and carry_free(multiplier*t)) for t in short)
        proofs[gap] = ('high', height, len(short))
    residual = {g for g, row in proofs.items() if row[0] == 'open'}
    assert residual == {1356, 1424, 2936, 3364}
    assert not any(a + b in residual for a in residual for b in residual)
    return proofs


def low_difference_count():
    power = 5**9
    period = power // 5
    _, prefixes = low_prefixes(1, 9)
    allowed = set(prefixes)
    x, logs = 2, []
    for j in range(period):
        if x in allowed:
            logs.append(j)
        x = 16*x % power
    assert x == 2 and len(logs) == len(allowed) == 1145
    differences = {(i-j) % period for i in logs for j in logs}
    assert len(differences) == 151047
    return len(differences), period


def angular_gap_bound():
    # Round the allowed ratio intervals outward on a fixed grid.
    height = 7
    scale = 5**(height-1)
    cells = [(6*scale, 30*scale)]
    for i in range(height):
        cells = [pair for a, b in cells for pair in clip_bands(a, b, scale//5**i, 0)]
    assert len(cells) == 636
    unit = 1 << 36
    intervals = []
    for a, b in cells:
        for c, d in cells:
            for lift in (1, 5):
                lo = max(unit, lift*a*unit // d)
                hi = min(5*unit, (lift*b*unit + c - 1) // c)
                if lo <= hi:
                    intervals.append((lo, hi))
    end = unit
    lower = 0
    for lo, hi in sorted(intervals):
        if end < lo:
            # log(lo/end) >= (lo-end)/lo.
            lower += (lo-end)*unit // lo
        end = max(end, hi)
    assert end == 5*unit
    # log(5) < 161/100. The omitted fraction exceeds 39/100.
    assert lower * 10000 > 39*161*unit
    return lower, unit


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('certificate', type=Path, help='output JSON for the independent checker')
    args = parser.parse_args()
    cert = gap_certificates()
    args.certificate.write_text(json.dumps(cert, indent=1) + '\n')
    for kind in ('low', 'high', 'open'):
        print(kind, sum(value[0] == kind for value in cert.values()))
    print('Allowed finite differences modulo 5^8:', low_difference_count())
    print('Lower bound for omitted logarithmic length:', angular_gap_bound())
