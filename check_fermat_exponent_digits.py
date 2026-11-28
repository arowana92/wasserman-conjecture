"""Exact path counts for the forced cut at 13 and 65537. Python 3.9+."""

from itertools import accumulate
from math import comb, gcd

P = 65537
X = 49153
Y = P - X


def predecessor_step(base, x):
    """Act on continuation counts indexed by the preceding digit plus one."""
    inverse = pow(x, -1, base)
    permutation = [inverse * s % base - 1 for s in range(1, base)]
    intervals = []
    for i in range(1, base):
        left = x - x * i // base - 1
        intervals.append((left, left + i))

    def step(values):
        partial = list(accumulate((values[k] for k in permutation), initial=0))
        return [partial[right] - partial[left] for left, right in intervals]

    return step


def weight_transitions():
    """Reachable ordinary multiplication carries, seeded by the first digit 3."""
    base = 13
    start = (3 * X // base, (3 * Y + 1) // base)
    states = [start]
    indices = {start: 0}
    rows = []
    for a, b in states:
        children = []
        for d in range(base):
            next_a, low_a = divmod(X * d + a, base)
            next_b, low_b = divmod(Y * d + b, base)
            if low_a + low_b >= base:
                continue
            child = (next_a, next_b)
            if child not in indices:
                indices[child] = len(states)
                states.append(child)
            children.append(indices[child])
        rows.append(children)
    return rows


def predecessor_counts():
    step = predecessor_step(P, X)
    v = [1] * (P - 1)
    counts = [1]
    for depth in range(1, 10):
        w = step(v)
        counts.append(w[4])
        if depth == 9:
            assert all(8 * new <= (8 * 32768 + 1) * old for old, new in zip(v, w))
        v = w
    assert counts[:4] == [1, 5, 196596, 6442549034]
    assert (8 * 32768 + 1) ** 16 < 8 ** 16 * P ** 15
    return counts


def weight_counts():
    rows = weight_transitions()
    assert len(rows) == 65536
    assert sum(map(len, rows)) == 458752
    v = [1] * len(rows)
    counts = [1]
    for depth in range(1, 21):
        w = [sum(v[j] for j in row) for row in rows]
        counts.append(w[0])
        if depth == 20:
            assert all(10**6 * new <= 7000001 * old for old, new in zip(v, w))
        v = w
    assert counts[:5] == [1, 6, 40, 286, 1956]
    assert counts[10] == 231470218
    assert 7000001**1000 < 13**759 * 10**6000
    return counts


def clear_layers(base, x, y, scale, layers):
    """Use integer quotients (Legendre defects), independently of the graph."""
    a = x * scale
    b = y * scale + 1
    modulus = 1
    for _ in range(layers):
        modulus *= base
        if a // modulus + b // modulus < (a + b) // modulus:
            return False
    return True


def direct_checks(counts):
    for depth in range(6):
        assert sum(clear_layers(13, X, Y, 3 + 13 * t, depth + 1)
                   for t in range(13**depth)) == counts[depth]
    # An accepted finite prefix in the predecessor base can be padded by zeros.
    for base in (5, 7, 8, 9, 11, 16):
        for x in range(2, base - 1):
            if gcd(x, base) != 1:
                continue
            low = pow(x, -1, base)
            step = predecessor_step(base, x)
            v = [1] * (base - 1)
            for depth in range(4):
                count = sum(clear_layers(base, x, base - x, low + base * t, depth + 2)
                            for t in range(base**depth))
                assert count == v[low]
                v = step(v)
    for w in range(3, 300, 13):
        # Scaled-down slopes retain the leading congruences at 13.
        assert clear_layers(13, 13, 4, w, 5) == (comb(17 * w + 1, 13 * w) % 13 != 0)


def synchronized_exponents():
    # 16 is 3 mod 13, of order 3 there. Removing its torsion part
    # gives the step for b increasing by four at a fixed low digit W=3.
    modulus = 13**21
    tau = pow(3, 13**20, modulus)
    step = 16 * pow(pow(16, 13**20, modulus), -1, modulus) % modulus
    w = tau
    last = []
    for b in range(4, 1000001, 4):
        w = w * step % modulus
        if b < 72:
            continue
        if clear_layers(13, X, Y, w, 20):
            last.append(b)
            assert not clear_layers(13, X, Y, w, 21)
    assert last == [50088, 248052]
    # Independently enumerate b/4 by residues and power 2 directly. Test
    # individual digit sums instead of quotient defects.
    roots = [1, tau, tau**2 % modulus]

    def passes(k, n):
        m = 13**n
        power = pow(16, k, m)
        r, = [root for root in roots if root * power % 13 == 3]
        z = r * power % m
        a, b = X * z, Y * z + 1
        for _ in range(n):
            a, left = divmod(a, 13)
            b, right = divmod(b, 13)
            if left + right >= 13:
                return False
        return True

    survivors = []
    for residue in range(13**3):
        if not passes(residue, 4):
            continue
        for k in range(residue, 250001, 13**3):
            if k >= 18 and passes(k, 20):
                survivors.append(4 * k)
                assert not passes(k, 21)
    assert sorted(survivors) == last
    assert pow(2, 12, 169) == 40 and (2**32 - 1) % P**2 == P * (P - 2)
    w = pow(3, 169, 2197)
    assert (w, X * w % 2197, (Y * w + 1) % 2197) == (1160, 936, 1391)


if __name__ == "__main__":
    c = predecessor_counts()
    d = weight_counts()
    direct_checks(d)
    synchronized_exponents()
    print("Predecessor counts through depth 3:", c[:4])
    print("Weight-prime counts at depths 1, 2, 4, 10:", [d[i] for i in (1, 2, 4, 10)])
    print("Joint count and period at depths (2, 10):", c[2] * d[10], 96 * P**2 * 13**10)
    print("All coordinate bounds, direct checks, and synchronized exponents passed.")
