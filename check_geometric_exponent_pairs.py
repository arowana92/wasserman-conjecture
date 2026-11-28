"""Independent counts of 60 low output bits and modular checks of exponent lifts."""

from collections import defaultdict
from math import gcd, isqrt

from carry_blocks import CarryGraph


def prime_power(n):
    for p in range(2, isqrt(n) + 1):
        if n % p == 0:
            while n % p == 0:
                n //= p
            return n == 1
    return n > 1


def multiplier_count(x, y, n):
    graph = CarryGraph(2, (x, y), (0, 0))
    # All endpoints work since only the bottom n columns are prescribed.
    words = [1] * len(graph.states)
    for _ in range(n - 1):
        words = [sum(words[t] for _, t in edges) for edges in graph.edges]
    return sum(words[t] for bit, t in graph.edges[0] if bit == 1)


def equation_count(x, y, n):
    even, odd = (x, y) if x % 2 == 0 else (y, x)
    # The output of the odd slope has low bit 1.
    counts = {-even//2: 1}
    for _ in range(1, n):
        following = defaultdict(int)
        for quotient, multiplicity in counts.items():
            for delta in (0, odd, -even):
                numerator = quotient + delta
                if numerator % 2 == 0:
                    following[numerator//2] += multiplicity
        counts = following
    return sum(counts.values())


def test_table():
    cuts = total = 0
    expected = {31: (92, 4127843005161), 100: (583, 19462275734923),
                250: (3218, 85421792224053)}
    for d in range(3, 251):
        if d % 2 == 1 and prime_power(d):
            for x in range(1, (d+1)//2):
                y = d - x
                if gcd(x, d) != 1 or (x & (x-1) == 0 and y & (y-1) == 0):
                    continue
                count = multiplier_count(x, y, 60)
                assert count == equation_count(x, y, 60), (d, x)
                if d <= 31:
                    mask = 2**10 - 1
                    direct = sum((x*t & y*t & mask) == 0 for t in range(1, mask+1, 2))
                    assert equation_count(x, y, 10) == direct
                cuts += 1
                total += count
        if d in expected:
            assert (cuts, total) == expected[d]
            print(f"D <= {d}: {cuts} cuts, at most {total} exceptional odd residues", flush=True)
    assert total*100000 < 15*2**59


def test_lifts():
    cases = 0
    for p in (2, 3, 5, 7):
        for q in (2, 3, 5, 7, 17, 37):
            if p == q:
                continue
            order = 1
            while pow(q, 2*order, p) != 1:
                order += 1
            c = q**order
            assert c % p in (1, p-1)
            if p > 2 and c % p != 1:
                for k in range(5):
                    layer = p**(k+1)
                    for e in (1, 3, 37):
                        t = pow(c, e*p**k, layer)
                        assert t == layer-1
                        for d in range(2, min(22, layer)+1):
                            assert all(x*t % layer + (d-x)*t % layer >= layer
                                       for x in range(1, d))
            a = c*c if p == 2 or c % p != 1 else c
            s = 0
            unit = a - 1
            while unit % p == 0:
                s += 1
                unit //= p
            for n in (2, 3, 4):
                modulus = p**n
                for offset in (1, 10**11 + 127):
                    previous = None
                    for k in range(n + 4):
                        low = p**(s+k)
                        values = [pow(a, e*p**k, low*modulus)
                                  for e in range(offset, offset + modulus)]
                        assert all(t % low == 1 for t in values)
                        residues = [t//low for t in values]
                        assert len(set(residues)) == modulus
                        if previous is not None:
                            level = p**min(n, s+k-2)
                            assert all((u-v) % level == 0 for u, v in zip(residues, previous))
                        previous = residues
                        cases += len(values)
    print(f"{cases} normalized exponent residues checked", flush=True)


if __name__ == "__main__":
    test_lifts()
    test_table()
