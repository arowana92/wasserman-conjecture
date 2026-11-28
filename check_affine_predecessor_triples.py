"""Exact controls for affine_predecessor_triples.md. Uses the standard library."""

from collections import Counter
from fractions import Fraction
from itertools import combinations, product, zip_longest
from math import gcd, isqrt, lcm, prod


def check_determinants():
    forms = list(product((1, 2, 3, 4), (-3, -1, 1, 2, 4)))
    counts = Counter()
    for triple in combinations(forms, 3):
        ordered = sorted(triple, key=lambda ab: Fraction(-ab[1], ab[0]))
        a, b = zip(*ordered)
        det = [[a[i] * b[j] - a[j] * b[i] for j in range(3)] for i in range(3)]
        if any(det[i][j] == 0 for i, j in combinations(range(3), 2)):
            continue
        bound = max(abs(n) for row in det for n in row)
        ds = [prod(d for d in row if d) for row in det]
        for t in range(21, 101):
            ws = [s * t + h for s, h in ordered]
            if any(gcd(x, y) > 1 for x, y in combinations(ws, 2)):
                continue
            m = prod(ws)
            ks = [pow(m // w, -1, w) for w in ws]
            g = gcd(*ks)
            counts["coprime triples"] += 1
            counts["gcd at least three"] += g >= 3
            row_sum = sum(k * (m // w) for k, w in zip(ks, ws))
            counts["row sum and gcd at least three"] += row_sum == m + 1 and g >= 3
            ms = []
            for i, w in enumerate(ws):
                mi, rem = divmod(ds[i] * ks[i] - a[i] ** 2, w)
                assert rem == 0
                assert (0 <= mi < ds[i]) if ds[i] > 0 else (ds[i] <= mi < 0)
                ms.append(mi)
            assert g <= bound**5 + 2 * max(a)**3 * bound**2
            if ms[0] == 0:
                assert a[0] ** 2 % g == 0
            else:
                resultant = (ms[0] * ms[1] * det[0][1] + ms[0] * a[0] * a[1]**2
                             - ms[1] * a[1] * a[0]**2)
                assert resultant != 0 and resultant % g == 0
    return dict(counts)


def value(poly, t):
    result = 0
    for coefficient in reversed(poly):
        result = result * t + coefficient
    return result


def v2(x):
    assert x != 0
    return (x & -x).bit_length() - 1


def check_prefix_counts():
    # Coefficients in increasing order. The second row has two quadratic
    # parts and a cubic total whose derivative has two integer roots.
    rows = [((0, 1), (0, 2), (0, 4)),
            ((0, 2, 3, 1), (2, 3, 1), (0, 4, 2)),
            ((0, 1, 0, 1), (7, 1, 2, 4), (2, 1)),
            ((6, 11, 6, 1), (0, 2, 1), (0, 3, 4, 1))]
    balls = prefixes = passing = 0
    # 14655/10000 < lambda, by substitution into X^3-X^2-1.
    assert 14655**3 < 14655**2 * 10000 + 10000**3
    for parts in rows:
        total = tuple(map(sum, zip_longest(*parts, fillvalue=0)))
        polys = (total,) + parts
        derivs = [tuple(j * c for j, c in enumerate(f) if j) for f in polys]
        for n in (5, 8, 11):
            pending = [(0, 0)]
            while pending:
                x, h = pending.pop()
                ds = [value(f, x) for f in derivs]
                if any(d % (1 << h) == 0 for d in ds):
                    if h < n:
                        pending.extend(((x, h + 1), (x + (1 << h), h + 1)))
                    continue
                es = [v2(d) for d in ds]
                count = 0
                for t in range(x, 1 << n, 1 << h):
                    # Exactly the low n+e_i bits are fixed by this input residue.
                    # Check collisions directly, without reproducing the branching proof.
                    used = 0
                    for i, (f, e) in enumerate(zip(polys, es)):
                        mask = (1 << (n + e)) - 1
                        v = value(f, t)
                        assert v2(value(f, t + (1 << n)) - v) == n + e
                        digits = ((~v if i == 0 else v) & mask) >> (h + min(es))
                        if used & digits:
                            break
                        used |= digits
                    else:
                        count += 1
                exponent = n - h + max(es) - min(es)
                assert count * 10000**exponent <= 14655**exponent
                balls += 1
                prefixes += 1 << (n - h)
                passing += count
    return balls, prefixes, passing


def is_prime(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))


def check_family():
    checked = 0
    for t in range(5, 1001):
        w, z = 5 * t + 3, 8 * t + 5
        for j in (0, 1, 2, 12345):
            y = 2 * w + z + j * w * z
            # Compute the idempotents directly, including at composite inputs.
            assert pow(z * y, -1, w) == 25
            assert pow(w * y, -1, z) == 32
            checked += 1
    w, z, y = 5007713, 8012341, 18027767
    assert all(is_prime(n) for n in (w, z, y))
    assert 5 * z - 8 * w == 1 and y == 2 * w + z
    assert 1 + w * z * y == 723337184691130355012
    period = lcm(w - 1, z - 1)
    assert period == 10030872791520 and pow(y, period, w * z) == 1
    for j in (0, 1, 2, 10**6):
        power = pow(y, 1 + period * j, w * z)
        assert (pow(z * power, -1, w), pow(w * power, -1, z)) == (25, 32)
    return checked, period


if __name__ == "__main__":
    print("CRT checks", check_determinants())
    print("Stopped balls, prefixes, collision-free prefixes", check_prefix_counts())
    print("Family checks and exponent period", check_family())
