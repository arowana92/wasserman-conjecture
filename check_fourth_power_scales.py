"""Unbounded cofactor-35 scale checks and the two surviving prime inputs."""

from collections import deque
from math import gcd, isqrt


def clear(parts, modulus):
    return sum(v % modulus for v in parts) < modulus


def ternary_parts(c):
    q = 3**(4*c)
    # Digits of r = s^4 satisfying the independent cut 16r + 19r.
    stack = [(1, 3)]
    while stack:
        r, m = stack.pop()
        if m == q:
            z = pow(35*r, -1, q)
            if z % 4 == 0 and 35*z < 19*q:
                yield z
            continue
        for digit in range(3):
            v = r + digit*m
            if clear((16*v, 19*v), 3*m):
                stack.append((v, 3*m))


def check_binary_shape(q, z):
    slopes = (16*q, 35*z, 19*q - 35*z)
    # For r = 1 modulo 16, start ordinary multiplication with its low four bits.
    seed = (slopes[0], slopes[1], slopes[2] + 1)
    if any(not clear(seed, m) for m in (2, 4, 8, 16)):
        return 0
    # Cheap finite prefix filter; dropping a branch requires an actual carry.
    prefixes = [1]
    m = 16
    for _ in range(8):
        prefixes = [v for u in prefixes for v in (u, u + m)
                    if clear((slopes[0]*v, slopes[1]*v, slopes[2]*v+1), 2*m)]
        m *= 2
        if not prefixes:
            return 0
    start = tuple(v // 16 for v in seed)
    queue = deque([start])
    seen = {start}
    while queue:
        a, b, c = queue.popleft()
        # If the unread digits are zero the multiplication carries are exactly
        # the remaining outputs. Any pairwise-disjoint triple would accept.
        assert (a & b) or (a & c) or (b & c)
        for x, y, z1 in ((a, b, c), (a+slopes[0], b+slopes[1], c+slopes[2])):
            if (x & 1) + (y & 1) + (z1 & 1) > 1:
                continue
            state = (x // 2, y // 2, z1 // 2)
            if state not in seen:
                seen.add(state)
                queue.append(state)
    return len(seen)


def check_small_exponents():
    table = []
    for c in range(1, 6):
        q = 3**(4*c)
        zs = list(ternary_parts(c))
        if c <= 3:
            # Reverse the enumeration. Read z first, then check the scale by
            # quotient defects at every power of 3 through q.
            cross = set()
            for z in range(4, (19*q-1)//35+1, 4):
                if gcd(z, 3) != 1:
                    continue
                r = pow(35*z, -1, q)
                if r % 3 == 1 and all(16*r//3**k + 19*r//3**k == 35*r//3**k
                                      for k in range(1, 4*c+1)):
                    cross.add(z)
            assert cross == set(zs)
        sizes = [check_binary_shape(q, z) for z in zs]
        table.append((len(zs), sum(v != 0 for v in sizes), max(sizes)))
    assert table == [(2,0,0),(26,2,19),(301,23,3131),(4501,413,7955),(71002,6748,10190)]
    print("Three-adic exponents 1 through 5:", table)


def prime_certificate(p, witness, factors):
    product = 1
    for ell, power in factors:
        assert all(ell % d for d in range(2, isqrt(ell)+1))
        assert pow(witness, (p-1)//ell, p) != 1
        product *= ell**power
    assert product == p-1 and pow(witness, p-1, p) == 1


def check_constants():
    prime_certificate(21523361, 6, [(2,5),(5,1),(17,1),(41,1),(193,1)])
    prime_certificate(926510094425921, 7, [(2,6),(5,1),(17,1),(41,1),(193,1),(21523361,1)])
    for n, divisor in ((8,17),(12,41),(20,41),(24,17),(28,41)):
        s = (3**n+1)//2
        assert s % divisor == 0 and gcd(s//divisor, divisor) == 1


def check_exponent_interval(s, limit):
    c0, = [c for c in range(1,7) if 3**c*s % 8 in (3,5) and 3**c*s % 7 in (3,4)]
    end = (limit-c0)//6
    step = 81**6
    mask = (1 << 110)-1
    initial = s**4 * 81**c0
    advances = [pow(step, 1 << j, mask+1) for j in range(end.bit_length()+1)]
    tables = []
    for p in (5, 7):
        period = p**6
        modulus = p*period
        w = initial % modulus
        allowed = bytearray(period)
        powers = [p**i for i in range(1,8)]
        for i in range(period):
            allowed[i] = all(clear((16*w,19*w+1), m) for m in powers)
            w = w*step % modulus
        assert w == initial % modulus
        tables.append(allowed)
    stack = [(0,initial & mask,0)]
    terminals = static = survivors = 0
    while stack:
        i, w, depth = stack.pop()
        if 16*w & (19*w+1) & ((1 << (depth+5))-1):
            continue
        if i+(1 << depth) > end:
            terminals += 1
            if all(a[i % len(a)] for a in tables):
                static += 1
                # Recompute terminal residues directly rather than trusting
                # the recursively advanced scale or the lookup tables.
                direct = s**4 * pow(81, c0+6*i, mask+1) & mask
                assert direct == w
                for p in (5,7):
                    z = s**4 * pow(81, c0+6*i, p**7)
                    assert all(16*z//p**j + (19*z+1)//p**j == (35*z+1)//p**j
                               for j in range(1,8))
                if not (16*w & (19*w+1) & mask):
                    survivors += 1
            continue
        stack.append((i+(1 << depth), w*advances[depth] & mask, depth+1))
        stack.append((i,w,depth+1))
    return terminals, static, survivors


if __name__ == "__main__":
    check_constants()
    check_small_exponents()
    for n, expected in ((16,(5383615,17548,0)), (32,(5384016,17647,0))):
        s = (3**n+1)//2
        result = check_exponent_interval(s, 10**12)
        assert result == expected
        print("Prime and exponent counts:", s, result)
