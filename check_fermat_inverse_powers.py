"""Independent finite checks for the two bounded inverse lifts at 65537."""

from fractions import Fraction
from collections import defaultdict
from itertools import combinations
from math import gcd, isqrt, prod

P = 65537


def factor(n):
    factors = {}
    for d in range(2, isqrt(n) + 1):
        while n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            n //= d
    if n > 1:
        factors[n] = 1
    return factors


def candidates(h, exponent):
    # First traverse the cyclic subgroup and use its order at each odd prime.
    order = (P - 1) // exponent
    generator = pow(3, exponent, P)
    group = {pow(generator, j, P) for j in range(order)}
    assert len(group) == order
    rows = []
    for a in sorted(group):
        for x in range(a, h * P, P):
            good = []
            for p in factor(x):
                if p > 2 and pow(-h * P, (p - 1) // gcd(p - 1, exponent), p) == 1:
                    good.append(p)
            if good:
                rows.append((x, tuple(good)))
    return sorted(rows)


def candidates_by_enumeration(h, exponent):
    # List actual power residues at P, then all potential inputs at each divisor.
    group = {pow(t, exponent, P) for t in range(1, P)}
    output = []
    for x in range(1, h * P):
        if x % P == 0 or pow(x, -1, P) not in group:
            continue
        primes = tuple(p for p in factor(x)
                       if p > 2 and any((h * P * pow(t, exponent, p) + 1) % p == 0
                                        for t in range(1, p)))
        if primes:
            output.append((x, primes))
    return output


def check_odd_companion():
    slopes = (16384, 32768, 16129)
    total = sum(slopes)
    live = [tuple(x // 2**j for x in slopes) for j in range(1, 17)]
    assert live[-1] == (0, 0, 0) and len(set(live)) == 16
    edges = set()
    for target in live:
        for digit in (0, 1):
            for units in ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)):
                src = tuple(2 * c + u - digit * x for c, u, x in zip(target, units, slopes))
                bounds = list(zip(src + (sum(src),), slopes + (total,)))
                low = max([Fraction(0)] + [Fraction(c, x) for c, x in bounds])
                high = min([Fraction(1)] + [Fraction(c + 1, x) for c, x in bounds])
                if low < high:
                    assert src in live
                    edges.add((live.index(src), digit, live.index(target)))
    expected = {(i, 0, min(i + 1, 15)) for i in range(16)} | {(15, 1, 0)}
    assert edges == expected

    # Independent forward search above the forced prefix 257 modulo 16384.
    initial = tuple(x * 257 // 16384 for x in slopes)
    stack = [initial]
    visited = {initial}
    while stack:
        a, b, c = state = stack.pop()
        assert (a & b) or (b & c) or (c & a)
        for digit in (0, 1):
            raw = tuple(s + digit * x for s, x in zip(state, slopes))
            if sum(v % 2 for v in raw) > 1:
                continue
            nxt = tuple(v // 2 for v in raw)
            if nxt not in visited:
                visited.add(nxt)
                stack.append(nxt)
    assert len(visited) == 3201
    assert pow(65281, -1, 16384) == 257
    return len(edges), len(visited)


def check_spacing():
    cases = sharp = 0
    for h in range(3, 21):
        a = 2**(h - 2)
        for b in range(1, h - 1):
            slopes = (a, 2 * a, a - 2**b + 1)
            if 2 * b > h:
                z = 1 + 2**(h - b)
                x, y, w = (s * z for s in slopes)
                assert ((x & y) | (x & w) | (y & w)) == 0
                sharp += 1
                continue
            for z in (*range(1, 1025), 1 + 2**h, 1 + 2**(h - 1), 1 + 2**h + 2**(2*h)):
                x, y, w = (s * z for s in slopes)
                clear = ((x & y) | (x & w) | (y & w)) == 0
                positions = [j for j in range(z.bit_length()) if z >> j & 1]
                assert clear == all(j - i >= h for i, j in zip(positions, positions[1:]))
                cases += 1
    return cases, sharp


def smaller_exponent_companions():
    u, a = 65281, 8192
    prefix = pow(u, -1, a)
    count = forward = reverse = 0
    for i, j in combinations(range(1, u // a + 1), 2):
        if (i + j) * a >= u:
            continue
        slopes = (i * a, j * a, u - (i + j) * a)
        start = tuple(x * prefix // a for x in slopes)
        todo = [start]
        seen = {start}
        while todo:
            s = todo.pop()
            x, y, z = s
            assert (x & y) or (y & z) or (x & z)
            for digit in (0, 1):
                raw = tuple(c + digit * f for c, f in zip(s, slopes))
                if sum(n % 2 for n in raw) > 1:
                    continue
                dest = tuple(n // 2 for n in raw)
                if dest not in seen:
                    seen.add(dest)
                    todo.append(dest)
        forward += len(seen)

        # Independently enumerate rational cells, including the weight floor.
        augmented = slopes + (u,)
        cells = {tuple(x * n // den for x in augmented)
                 for den in augmented for n in range(den)}
        clear = {s for s in cells if sum(s[:3]) == s[3]}
        parents = defaultdict(list)
        for s in clear:
            for digit in (0, 1):
                dest = tuple((v + digit * x) // 2 for v, x in zip(s, augmented))
                if dest in clear:
                    parents[dest].append(s)
        zero = (0, 0, 0, 0)
        stack = [zero]
        live = {zero}
        while stack:
            for src in parents[stack.pop()]:
                if src not in live:
                    live.add(src)
                    stack.append(src)
        assert start + (u * prefix // a,) not in live
        reverse += len(cells)
        count += 1
    assert count == 9
    return count, forward, reverse


def check():
    assert factor(P) == {P: 1}
    for h, exponent, target in [
        (1, 16384, [(65281, (97, 673))]),
        (1, 8192, [(65281, (97, 673))]),
        (2, 8192, [(65521, (65521,)), (65793, (7,))]),
    ]:
        assert candidates(h, exponent) == candidates_by_enumeration(h, exponent) == target
        m = 4 * exponent
        for x, primes in target:
            assert x % 2 and (x % m) & ~(h * P % m)  # even t forces a binary carry
            g = prod(p**factor(x)[p] for p in primes)
            if h == 2:
                assert (x % m) & ~((2 * P + 1) % m)
                assert (2 * P + 1) % (2 * m) == 3
                assert g <= 2 * m
            elif exponent == 16384:
                assert (P + 1) % m == 2 and g <= m
            else:
                assert g == 65281 and all(d <= m for d in (1, 97, 673))
        print(h, exponent, target)
    # Prime classes giving independent, unrestricted exponents.
    first, second = 2359333, 23855467
    moduli = (3, 11, 331, 43691)
    assert all(factor(p) == {p: 1} for p in (first, second, *moduli))
    assert first % P == 1 and second % P == P - 1
    assert all(pow(q, (p - 1) // 2, p) == 1 for q in (first, second) for p in moduli)
    assert factor(P + 1) == {2: 1, 3: 2, 11: 1, 331: 1}
    assert factor(2 * P - 1) == {3: 1, 43691: 1}
    assert all(pow(2, (p - 1) // 2, p) == p - 1 for p in (3, 11, 331))
    assert all(p % 4 == 3 for p in (3, 43691))
    candidates4096 = [(49153, (13,)), (65281, (97, 673)), (65533, (13,))]
    assert candidates(1, 4096) == candidates_by_enumeration(1, 4096) == candidates4096
    assert [x for x, _ in candidates4096 if (x % 16384) & ~(P % 16384) == 0] == [49153]
    assert 49153 * 4 % P == 1 and (3 * P + 1) % 13 == 0
    assert (49153 % 32768) & ~(P % 32768)


if __name__ == "__main__":
    check()
    print("Companion graph transitions and forward states:", check_odd_companion())
    print("Spacing equivalences and sharpness witnesses:", check_spacing())
    print("4096th-power companion directions, forward states, reverse cells:",
          smaller_exponent_companions())
