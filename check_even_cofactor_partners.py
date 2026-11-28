"""Verify the even two-prime cofactors through 32768 and the three remaining rows."""

from collections import Counter, defaultdict
from fractions import Fraction
from math import gcd, prod

from carry_blocks import accepted_hull, intersection
from check_even_two_power_cofactors import (merge_progressions, positions, prefix_graph,
                                           primes_and_powers, verify_shape)
from check_even_two_power_cross import reconstruct


# w, P, u, a, b, initial residue, refined residue, refined modulus
CASES = ((64, 241, 3, 51, 6, 24491, 579755, 89952768),
         (32, 761, 5, 15, 10, 33767, 3483219159047, 5106984550400),
         (256, 109, 11, 143, 33, 92107, 73612323019, 160927055872))


def extension():
    primes, powers = primes_and_powers(32768 // 32)
    counts = Counter()
    remaining = []
    for j in range(5, 14):
        w = 2**j
        selected = [(P, p) for P, p in powers if 15423 < w*P <= 32768]
        actual = defaultdict(set)
        for P, p, u, a, b, r, modulus, slopes in positions(w, selected, primes, 32768):
            actual[P].add((w, P, u, a, b))
            reason = verify_shape(w, P, p, u, a, b, r, modulus, slopes)
            counts[reason] += 1
            if reason is None:
                remaining.append((w, P, u, a, b, r))
        for P, _ in selected:
            assert actual[P] == reconstruct(w * P)
    assert sorted(remaining) == sorted(case[:6] for case in CASES)
    assert sum(counts.values()) == 217370
    print("Additional assignments", dict(counts), flush=True)


class FixedPrefix:
    """Carry language with a prescribed low block and only states leading to acceptance."""

    def __init__(self, base, power, residue, suffix):
        self.base = base
        self.power = power
        self.low = residue % power
        self.suffix = suffix
        parents = [[] for _ in suffix.edges]
        for i, row in enumerate(suffix.edges):
            for _, j in row:
                parents[j].append(i)
        live = {i for i in range(len(parents)) if suffix.accept >> i & 1}
        todo = list(live)
        while todo:
            for i in parents[todo.pop()]:
                if i not in live:
                    live.add(i)
                    todo.append(i)
        assert 0 in live
        self.edges = []
        low = self.low
        while power > 1:
            self.edges.append([(low % base, len(self.edges) + 1)])
            low //= base
            power //= base
        offset = len(self.edges)
        ids = {i: j + offset for j, i in enumerate(sorted(live))}
        self.edges.extend([(digit, ids[j]) for digit, j in suffix.edges[i] if j in ids]
                          for i in sorted(live))
        self.accept = sum(1 << ids[i] for i in live if suffix.accept >> i & 1)

    def accepts_integer(self, n):
        return n % self.power == self.low and self.suffix.accepts_integer(n // self.power)


def check_scales(case, digits):
    w, P, u, a, b, r, residue, modulus = case
    D = w*P
    slopes = (a*P, b*w, D - a*P - b*w)
    initial = (r, u*D)
    congruence = initial
    languages = []
    for base, layer in ((2, w), (u, u)):
        graph = prefix_graph(base, slopes, layer, r)
        h, step = accepted_hull(graph)
        congruence = merge_progressions(*congruence, r % layer + layer*h, layer*step)
        languages.append(FixedPrefix(base, layer, r, graph))
    assert congruence == (residue, modulus)
    # The last test deliberately uses only the weaker initial congruence.
    kept, nodes, direct = intersection(languages, 10**digits, initial[1], initial[0])
    assert not kept
    print(D, "refined class", congruence, "prefix intervals", nodes, "direct", direct,
          flush=True)


def factors(n):
    result = {}
    p = 2
    while p*p <= n:
        while n % p == 0:
            result[p] = result.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        result[n] = 1
    return result


def order(x, modulus, size):
    result = size
    for p in factors(size):
        while result % p == 0 and pow(x, result // p, modulus) == 1:
            result //= p
    assert pow(x, result, modulus) == 1
    return result


def roots_of_target(sizes, orders):
    density = Fraction(1)
    forbidden = []
    for ell in sorted(factors(prod(sizes))):
        ns = [factors(x).get(ell, 0) for x in sizes]
        bs = [factors(x).get(ell, 0) for x in orders]
        b = max(bs)
        if not b:
            continue
        h = min(n - v for n, v in zip(ns, bs) if v)
        forbidden.append(ell**(h+1))
        lifts = sum(ell**sum(min(n, j) for n in ns) for j in range(h+1))
        density *= Fraction((ell-1)*ell**(b-1)*lifts, ell**sum(ns))
    return density, forbidden


def prime_classes():
    individual = []
    for *_, r, m in CASES:
        sizes, orders = [], []
        for p, k in sorted(factors(m).items()):
            w = p**k
            if p == 2:
                assert r % 4 == 3
                sizes.extend((2, w//4))
                orders.extend((2, order(-r, w, w//4)))
            else:
                size = w//p*(p-1)
                sizes.append(size)
                orders.append(order(r, w, size))
        individual.append(roots_of_target(sizes, orders))
    a, b, c = (x for x, _ in individual)
    assert individual == [(Fraction(1,720), [2,3,5]), (Fraction(9,15200), [2,5,19]),
                          (Fraction(1,45), [2,9,5])]
    assert 1 - a - b - c + 4*a*c == Fraction(1201553,1231200)
    print("Possible prime fractions and excluded exponent divisors", individual)


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--digits", type=int, default=10000)
    parser.add_argument("--skip-extension", action="store_true")
    args = parser.parse_args()
    if not args.skip_extension:
        extension()
    for case in CASES:
        check_scales(case, args.digits)
    prime_classes()
