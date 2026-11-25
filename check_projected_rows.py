"""Check the determinant inequalities, binary encoding and prefix-seven offsets."""

from collections import Counter
from itertools import combinations, product
from math import gcd

from check_endpoint_pair_bound import prime_powers, rows

ORIGINS = ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1))
EDGES = tuple(combinations(range(4), 2))
MATCHINGS = (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2)))
DIRECTIONS = ((1, 1, 1), (1, 2, 4), (2, 3, 5))


def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


def geometry(row, z):
    shifts = [tuple(x-y for x, y in zip(row, e)) for e in ORIGINS]
    dets = {(i, j): abs(sum(x*y for x, y in zip(z, cross(shifts[i], shifts[j]))))
            for i, j in EDGES}
    eligible = all(any(cross(v, z)) for v in shifts)
    matches = [dets[e]*dets[f] for e, f in MATCHINGS if dets[e] and dets[f]]
    assert eligible == bool(matches)
    if eligible:
        zeros = [e for e in EDGES if not dets[e]]
        assert all(set(e).isdisjoint(f) for e, f in combinations(zeros, 2))
    return eligible, max(dets.values()), all(dets[0, i] for i in (1, 2, 3)), dets


def allocations(n, atoms, row, first, second):
    scales = [1]*4
    partial = False
    for p, w in atoms:
        left, right = gcd(w, first), gcd(w, second)
        if left > right:
            continue
        power = w//left
        partial |= 1 < power < w
        if power == 1:
            continue
        residue = tuple(x % power for x in row)
        i = ORIGINS.index(residue)
        assert (i == 0) == (n % p == 0)
        scales[i] *= power
    return scales, partial


def check_exact(limit):
    counts = Counter()
    for n in range(3, limit+1):
        m = n-1
        atoms = tuple(prime_powers(n))+tuple(prime_powers(m))
        table = [(r, k, [geometry(r, z) for z in DIRECTIONS]) for r, k in rows(n)]
        for r, first, data in table:
            for _, second, other in table:
                G = gcd(n*m, gcd(first, second))
                rs, partial = allocations(n, atoms, r, first, second)
                counts["ordered row pairs"] += 1
                counts["partially retained on first row"] += partial
                for eligible, height, stronger, dets in data:
                    assert all(dets[i, j] % (rs[i]*rs[j]) == 0 for i, j in EDGES)
                    if not eligible:
                        counts["proportional cases"] += 1
                        continue
                    counts["single-row bounds"] += 1
                    assert G**3*height**5 >= 27*m*m
                    counts["single-row bounds forcing a divisor"] += height**5 < 27*m*m
                    if stronger:
                        counts["stronger bounds"] += 1
                        assert G**6*height**9 >= 729*m**4
                    for valid, other_height, _, _ in other:
                        if valid:
                            counts["two-row bounds"] += 1
                            assert G*height**2*other_height**2 >= n*m
    return counts


def codeword(m, u, v, shift, length):
    x, y = u*m, v*m
    coded = 0
    for i in range(length):
        if (y << shift) & (1 << i):
            assert x & (1 << i) == 0
        else:
            coded |= y & (1 << i)
    return coded


def decode(coded, u, v, shift, length):
    m = 0
    for i in range(length):
        bit = 1 << i
        if ((v*m) << shift) & bit:
            multiplier, required = u, 0
        else:
            multiplier, required = v, coded & bit
        if (multiplier*m) & bit != required:
            m += bit
    return m


def check_encoding():
    length = 11
    mask = (1 << length)-1
    total = 0
    for u, v, shift in product(range(1, 16, 2), range(1, 16, 2), range(1, 5)):
        words = set()
        permitted = {c for c in range(mask+1) if c & (c << shift) == 0}
        for m in range(mask+1):
            if u*m & ((v*m) << shift) & mask:
                continue
            word = codeword(m, u, v, shift, length)
            assert word not in words and word in permitted
            assert decode(word, u, v, shift, length) == m
            words.add(word)
        if u == v == 1:
            assert words == permitted
        total += len(words)
    return total


def direct_seven(exponent, bits):
    # Allocate each low bit to a single part or leave it out of c altogether.
    q = 2**exponent
    assert bits <= exponent and (7*((1 << bits)-1))**5 < 27*(7*q)**2
    found = set()
    for allocation in product(range(4), repeat=bits):
        low = [0]*3
        for i, dest in enumerate(allocation):
            if dest:
                low[dest-1] += 1 << i
        row = tuple(a*q+b for a, b in zip((1, 2, 4), low))
        if gcd(*row) < 3:
            continue
        for e in ORIGINS:
            x, y, z = (a-b for a, b in zip(row, e))
            if y == 2*x and z == 4*x:
                found.add(sum(low))
                break
    return found


def spaced(m):
    return (m & ((m << 1) | (m << 2))) == 0


def spaced_strings(cap, prefix=0, bit=1):
    if bit > cap:
        yield prefix
        return
    yield from spaced_strings(cap, prefix, 2*bit)
    if prefix+bit <= cap:
        yield from spaced_strings(cap, prefix+bit, 8*bit)


def check_seven():
    small = direct_seven(32, 9)
    assert small == {7*m for m in range(512//7+1) if 7*m < 512 and spaced(m)}
    H = 10**6
    cap = H//7
    direct = {m for m in range(1, cap+1) if spaced(m)}
    generated = list(spaced_strings(cap))
    assert len(generated) == len(direct)+1 and set(generated) == direct | {0}
    k = 1
    while (7*H)**5 >= 27*(7*(1 << k))**2 or H >= (1 << k):
        k += 1
    print(f"Prefix 7: {len(small)-1} remaining positive offsets through 511; "
          f"{len(direct)} through {H}; uniform sufficient exponent {k}.")


if __name__ == "__main__":
    print(f"Exact coefficients: {dict(check_exact(36))}")
    print(f"Checked {check_encoding()} binary encodings and inverse images.")
    check_seven()
