"""Verify the prescribed-prime examples and the finite binary residue covering."""

from math import comb, gcd, lcm, prod

from check_primitive_row_scales import direct_layers, powers


def is_prime(p):
    return p >= 2 and powers(p) == [(p, p)]


def crt(congruences):
    r, m = 0, 1
    for t, n in congruences:
        if n == 1:
            continue
        r += m*(((t-r)*pow(m, -1, n)) % n)
        m *= n
    assert all(r % n == t % n for t, n in congruences)
    return r, m


def fill(modulus, gap, length):
    """Residues of bounded words, with the first gap entries representing length <= 0."""
    mask = (1 << modulus)-1
    layers = [1]*gap
    for j in range(length):
        shift = pow(2, j, modulus)
        base = layers[-gap]
        moved = ((base << shift) & mask) | (base >> (modulus-shift))
        layers.append(layers[-1] | moved)
    return layers


def word_for(layers, modulus, gap, target):
    j = len(layers)-gap
    assert (layers[-1] >> target) & 1
    word = 0
    while j > 0:
        if (layers[j+gap-2] >> target) & 1:
            j -= 1
            continue
        word |= 1 << (j-1)
        target = (target-pow(2, j-1, modulus)) % modulus
        j -= gap
    assert target == 0
    return word


def exact_example():
    n = 245623
    a = (35089, 70178, 140356)
    c = (25662, 109746, 110215)
    assert sum(a) == sum(c) == n and direct_layers(n, (a, c))
    ka, kc = (comb(n, x)*comb(n-x, y) for x, y, _ in (a, c))
    common = gcd(ka, kc, n*(n-1))
    assert ka % 2 == 1 and powers(common) == [(p, p) for p in (3, 7, 13, 47, 67)]
    return common


def concrete(v):
    n = 7*v
    assert v % 455175 == 234928 and v % 8 == 1 and not (v & (v << 1 | v << 2))
    r, rem = divmod(n-1, 510)
    assert rem == 0 and r % 255 == 37 and r % 7 == 1 and r % 2
    a = (v, 2*v, 4*v)
    c = (238*r, 210*r, 62*r+1)
    assert sum(a) == sum(c) == n and gcd(v, 7) == 1
    assert all(sum(z % w for z in c) == 1 for w in (15, 17, r))
    assert sum(z//7 % 7 for z in c) >= 7 and all(z % 7 == 0 for z in c)


def concrete_controls():
    L = 455175
    layers = fill(L, 3, 41)
    totals = [bin(s).count("1") for s in layers[-2:]]
    # Independent enumeration by the position of the final one.
    seen40 = bytearray(L)
    seen41 = bytearray(L)

    def visit(value, start):
        if value < 2**40:
            seen40[value % L] = 1
        seen41[value % L] = 1
        for i in range(start, 41):
            visit(value+2**i, i+3)

    visit(0, 0)
    assert totals == [sum(seen40), sum(seen41)] == [455152, 455175]
    checks = 0
    for k in (0, 1, 9, 25, 70):
        for w in (0, (2**k-1)//7):
            if w & (w << 1 | w << 2):
                continue
            low = 1+8*w+2**(k+48)
            target = ((234928-low)*pow(2**(k+5), -1, L)) % L
            z = word_for(layers, L, 3, target)
            v = low+(z << (k+5))
            assert 2**(k+48) <= v < 2**(k+49)
            concrete(v)
            checks += 1
    for k in (0, 1, 5):
        concrete((8**(99+4760*k)-1)//7)
        checks += 1
    return totals, checks


def general_controls(target, cutoff):
    primes = [p for p in range(3, cutoff+1, 2) if is_prime(p)]
    if target == 1:
        auxiliary = next(p for p in range(cutoff+1, 2*cutoff+1) if is_prime(p))
        target = auxiliary
        primes = [p for p in range(3, auxiliary+1, 2) if is_prime(p)]
    h = (target & -target).bit_length()-1
    a = 1
    while 3*a <= max(cutoff, primes[-1]):
        a *= 3
    if h:
        while a < max(cutoff, primes[-1], h+1):
            a *= 3
    F = 2**max(h, 1)
    z = 2**a
    u = 1+z+z*z
    assert gcd(u, prod(primes)) == 1
    # Use powers of two distinct auxiliary primes for d,e so that a common
    # multiple of the modular orders is easy to certify at large cutoffs.
    candidates = (p for p in range(primes[-1]+1, 10*primes[-1]+100)
                  if is_prime(p) and u % p)
    q1, q2 = next(candidates), next(candidates)

    def big_power(q):
        w = q
        while w <= u:
            w *= q
        return w

    d, e = big_power(q1), big_power(q2)
    D = d*e
    r0, _ = crt(((pow(u*e, -1, d), d), (pow(u*d, -1, e), e)))
    tests = [(1, u), ((1+F*D*r0)*pow(u, -1, D**2), D**2)]
    orders = [3*a, (d*d//q1)*(q1-1), (e*e//q2)*(q2-1)]
    included = []
    for p in primes:
        j, w = 0, 1
        while target % (w*p) == 0:
            w *= p
            j += 1
        if not j:
            tests.append((-pow(u, -1, p), p))
            orders.append(p-1)
            continue
        first = p
        while first <= u:
            first *= p
        modulus = first*w*w
        tests.append((-w, modulus))
        orders.append(modulus//p*(p-1))
        included.append((p, j, first))
    R, L = crt(tests)
    stride = 1+(2*a+1)*lcm(*orders)
    assert pow(2, stride, L) == 2
    # The following exponents encode an actual odd multiplier without
    # materializing it. No two differ by a or 2a.
    fixed = {0}
    if h:
        fixed.add(h)
        for i in range(h):
            start = h+2*a+1+i*(3*a+1)
            fixed.update((start, start+a))
    prefix = sum(1 << i for i in fixed)
    offset = max(fixed)+2*a+1
    tail = ((R-prefix)*pow(2**offset, -1, L)) % L
    bits = fixed | {offset+i*stride for i in range(tail.bit_length()) if tail >> i & 1}
    assert sum(pow(2, i, L) for i in bits) % L == R
    overlaps = (bits & {i+a for i in bits}) | (bits & {i+2*a for i in bits})
    assert len(overlaps) >= h if h else not overlaps
    # A small CRT representative agrees with the sparse number at all the
    # necessary odd moduli and through the prescribed low binary layer.
    V, _ = crt(((R, L), (prefix % (2**(h+1)), 2**(h+1))))
    n = u*V
    r, rem = divmod(n-1, F*D)
    c = (u*e*r, u*d*r, (F*D-u*(d+e))*r+1)
    if h:
        assert (n-1) % (2*F) == F and all(y%2 for y in c)
        assert sorted(v % F for v in (V, z*V, z*z*V)) == [0, 0, 1]
    assert rem == 0 and gcd(r, D) == 1 and sum(c) == n and min(c) > 0
    assert all(y % u == 0 for y in c) and gcd(u, V) == 1
    assert all(sum(y % w for y in c) == 1 for w in (d, e, r) if w > 1)
    assert 2**h*prod(p**j for p, j, _ in included) == target
    for p, j, first in included:
        w = p**j
        assert n % w == 0 and n % (w*p) and c[0] % p
        for h in range(j):
            layer = first*p**h
            residues = [(V//w)*t % layer for t in (1, z, z*z)]
            assert sum(residues) == 3*layer-u > layer
    assert all(n*(n-1) % p for p in primes if target % p)


if __name__ == "__main__":
    print("Exact endpoint gcd at 245623:", exact_example())
    print("Residue cover at 40 and 41 bits, constructed multipliers:", concrete_controls())
    tests = [(m, B) for m, minimum in ((1, 3), (3, 3), (5, 5), (7, 7), (9, 3),
                                        (15, 5), (2, 3), (64, 7), (360, 5), (3**7*5**4*19, 19), (97, 97))
             for B in (minimum, max(41, minimum))]
    for target, B in tests:
        general_controls(target, B)
    print("Prescribed small-prime parts:", len(tests))
