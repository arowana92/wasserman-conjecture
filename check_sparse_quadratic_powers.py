"""Three-part carry certificates at N=P(3P+2), with P a prime power."""
from itertools import combinations
from math import comb

from search_endpoint import atoms, smallest_factors


DYADIC_LIMIT = 1_000_000
KNOWN_SCALE = 200_000_000
RESIDUAL_ROWS = [(616,77,128), (872,32,109), (1432,179,288),
                 (1528,191,288), (1640,32,205), (1672,32,209),
                 (1808,113,320), (2072,64,259), (2168,271,448),
                 (2224,128,417), (2224,192,417), (2344,96,293)]


def possible_leading_digits(slopes):
    """Twelve carry-free binary columns inside the rational direction bounds."""
    total = sum(slopes)
    frontier = {(0, 0, 0)}
    for depth in range(1, 13):
        following = set()
        words = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
        if depth > 1:
            words.append((0, 0, 0))
        for old in frontier:
            for bits in words:
                weight = sum(bits)
                new = tuple(2*z + total*b - s*weight
                            for z, b, s in zip(old, bits, slopes))
                if all(s-total <= z <= s for s, z in zip(slopes, new)):
                    following.add(new)
        frontier = following
    return bool(frontier)


def excluding_depth(d, i, j, r):
    """Return a depth with no residue of Q that avoids every carry so far."""
    residues = [0]
    place = 1
    for depth in range(1, 201):
        if len(residues)*r > 2048:
            return None
        modulus = place*r
        third = pow(3, -1, modulus)
        following = []
        for start in residues:
            for x in range(start, modulus, place):
                if x % r == 0:
                    continue
                v = d*x+3
                total = (d*x+1)*v*third % modulus
                if v*x*i % modulus + v*x*j % modulus <= total:
                    following.append(x)
        if not following:
            return depth
        residues, place = following, modulus
    return None


def verify_low_cut(d, i, j, r, depth):
    """Check the rejected tree by factorial floors on actual positive rows."""
    stack = [(0, 1, 0)]
    while stack:
        prior, place, level = stack.pop()
        modulus = place*r
        for x in range(prior, modulus, place):
            if x % r == 0:
                continue
            q = x + modulus*(((-1-d*x)*pow(d*modulus, -1, 3)) % 3)
            p = (d*q+1)//3
            v = d*q+3
            parts = (v*q*i, v*q*j, v*(p-q*(i+j)))
            assert min(parts) > 0 and sum(parts) == p*v
            if p*v//modulus == sum(a//modulus for a in parts):
                assert level+1 < depth, (d, i, j, r, depth, q)
                stack.append((x, modulus, level+1))


def odd_power_cofactors():
    sieve = smallest_factors(2400)
    total = leading = cuts = 0
    survivors = []
    for d in range(2, 2401, 2):
        if d % 3 == 0:
            continue
        bases = [r for r, _ in atoms(d, sieve)]
        for i in range(1, d//6+1):
            for j in range(i, d//3-i+1):
                if i*j % d:
                    continue
                total += 1
                if (i & -i) == (j & -j):
                    continue
                if not possible_leading_digits((3*i, 3*j, d-3*(i+j))):
                    continue
                leading += 1
                for r in bases:
                    depth = excluding_depth(d, i, j, r)
                    if depth is not None:
                        verify_low_cut(d, i, j, r, depth)
                        cuts += 1
                        break
                else:
                    survivors.append((d, i, j))
    assert (total, leading, cuts) == (62935, 24302, 24290)
    assert survivors == RESIDUAL_ROWS
    for d, i, j in survivors:
        assert d % 8 == 0
        modulus = {1432: 179, 2072: 7}.get(d, 8)
        for root in range(1, modulus):
            if modulus == 8 and root % 2 == 0:
                continue
            q = root*root
            v = d*q+3
            residue = (d*q+1)*v*pow(3, -1, modulus) % modulus
            assert v*q*i % modulus + v*q*j % modulus > residue
    print("Through cofactor 2400:", total, "pairs,", leading, "past binary tests,",
          cuts, "verified residue exclusions. Residual shapes:", survivors)


def large_binary_blocks():
    singles = 0
    pairs = []
    for d in range(5, DYADIC_LIMIT+1, 2):
        if d % 3 == 0:
            continue
        high = (d*d-1)//3
        width = d.bit_length()
        middle = 4*d*pow(3, -1, 1 << width) % (1 << width)
        # Construct the possible i from their low bits. Higher bits cannot
        # change columns already checked in the products d*i and 3*i.
        candidates = [0]
        for bit in range(width):
            place = 1 << bit
            following = []
            for prior in candidates:
                for i in (prior, prior+place):
                    if 3*i < d and not ((d*i & ~high | 3*i & ~middle) & place):
                        following.append(i)
            candidates = following
            if candidates == [0] and 3*(2*place) >= d:
                break
        entries = [i for i in candidates if i and d*i & ~high == 0]
        singles += len(entries)
        for i, j in combinations(sorted(entries), 2):
            if 3*(i+j) < d and not (d*i & d*j or 3*i & 3*j):
                pairs.append((d, i, j))
    assert singles == 15433 and len(pairs) == 3
    assert all(i*j % d for d, i, j in pairs)
    print("Separated binary blocks:", singles, "single entries. Pairs:", pairs)


def bounded_binary_scales():
    sieve = smallest_factors(DYADIC_LIMIT)
    tested = 0
    for d in range(1, DYADIC_LIMIT+1, 2):
        if d % 3 == 0 or 2*d*d+1 <= 3*KNOWN_SCALE:
            continue
        # CRT idempotents select the powers at which 3*Q*i is 1.
        roots = [0]
        for _, power in atoms(d, sieve):
            f = d//power
            bit = f*pow(f, -1, power)
            roots += [(x+bit) % d for x in roots]
        q = 1 << ((2*d).bit_length()-1)
        while d*q+1 > 3*KNOWN_SCALE:
            if (d*q+1) % 3 == 0:
                tested += 1
                v = d*q+3
                target = d*((d*q+4)//3)
                inverse = pow(3*q, -1, d)
                for root in roots[1:]:
                    i = root*inverse % d
                    assert i and (i*(3*q*i-1)) % d == 0
                    assert 3*i >= d or v*i & ~target
            q //= 2
    assert tested == 1471770
    print("Overlapping binary blocks:", tested, "scales, zero single entries.")


def layer_controls():
    # Compare the middle-layer formula to the factorial floors, retaining
    # examples with no carry and with a carry in the tested layer.
    passed = failed = 0
    sieve = smallest_factors(10001)
    for p in range(5, 10001, 2):
        parts = (p*p, p*(p+1), p*(p+1))
        for r, w in atoms(p+1, sieve):
            if r == 2:
                continue
            b = (p+1)//w
            layer = r
            while layer <= w:
                modulus = w*layer
                free = sum(parts)//modulus == sum(a//modulus for a in parts)
                assert free == (4*(b % layer) > 3*layer)
                passed += free
                failed += not free
                layer *= r
    # Independently evaluate integrality and full parity for small rows.
    exact = 0
    for p in range(5, 65, 2):
        for q, w in atoms(3*p-1, sieve):
            if q != 2:
                continue
            v = 3*p+2
            for i in range(1, (p-1)//(2*w)+1):
                for j in range(i, (p-1)//w-i+1):
                    a, b = v*w*i, v*w*j
                    n = p*v
                    odd = comb(n, a)*comb(n-a, b) % 2
                    t = (n-1)//w
                    assert odd == (not (v*i & v*j or (v*i | v*j) & ~t))
                    exact += 1
    print("Middle-layer controls:", passed, "passing,", failed, "failing.",
          "Exact trinomial parity controls:", exact)


if __name__ == "__main__":
    layer_controls()
    odd_power_cofactors()
    large_binary_blocks()
    bounded_binary_scales()
