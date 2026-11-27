"""Decide the local bounded-digit endpoint tests by looking for infinite carry-free paths."""

from argparse import ArgumentParser
from fractions import Fraction
from functools import lru_cache

from carry_blocks import CarryGraph
from check_prime_power_cofactors import shapes
from search_endpoint import atoms, avoids, smallest_factors


def valuation(n, r):
    k = 0
    while n % r == 0:
        n //= r
        k += 1
    return k


@lru_cache(maxsize=128)
def infinite_graph(r, slopes, offsets):
    graph = CarryGraph(r, slopes, offsets)
    assert len(graph.states) <= 2*(sum(slopes)+1)
    # Starting with all states, retain those admitting paths one digit longer.
    live = set(range(len(graph.states)))
    while True:
        following = {i for i, edges in enumerate(graph.edges) if any(j in live for _, j in edges)}
        if following == live:
            break
        live = following
    return [{digit: j for digit, j in edges if j in live} for edges in graph.edges], live


def passes_prefix(r, slopes, offsets, x, k):
    edges, live = infinite_graph(r, slopes, offsets)
    state = 0
    for _ in range(k):
        x, digit = divmod(x, r)
        state = edges[state].get(digit)
        if state is None:
            return False
    return state in live


def required_rows(d, u, r, scalar):
    n = 1+d*scalar
    if u % r == 0:
        return ("A",) if valuation(n, r) == valuation(u, r) else ()
    if d % r == 0:
        return ("A", "C")
    if n % r == 0:
        return ("C",)
    return ("free",)


def locally_possible(d, shape, e, r):
    u, a, b, c, f, g, h = shape
    power_depth = 1+valuation(e, r)+(r == 2 and e % 2 == 0)
    k = max(power_depth, 1+valuation(u, r))
    modulus = r**k
    # Large primes outside the support of d U (d+1) accept scalar 1 outright.
    if d % r and u % r and (d+1) % r:
        return True
    scale_depth = k+valuation(d, r)
    scale_modulus = r**scale_depth
    for scalar in sorted({pow(t, e, modulus) for t in range(1, modulus) if t % r}):
        rows = required_rows(d, u, r, scalar)
        if "free" in rows:
            return True
        if "A" in rows and passes_prefix(r, (a, b, c), (0, 0, 1), scalar, k):
            return True
        if "C" in rows:
            v = (1+d*scalar)*pow(u, -1, scale_modulus) % scale_modulus
            if passes_prefix(r, (f, g, h), (0, 0, 0), v, scale_depth):
                return True
    return False


def find_first(e, limit):
    spf = smallest_factors(limit+1)
    count = 0
    for d in range(3, limit+1):
        factors = atoms(d, spf)
        for shape in shapes(d, factors):
            count += 1
            primes = {r for r, _ in factors+atoms(shape[0], spf)+atoms(d+1, spf)}
            if all(locally_possible(d, shape, e, r) for r in sorted(primes)):
                return d, shape, count
    return None, None, count


def check_247():
    shape = (3, 114, 117, 16, 1, 1, 1)
    for e, r in [(2, 3), (3, 13), (13, 13), (361, 19)]:
        assert not locally_possible(247, shape, e, r)
    assert {pow(t, 3, 13) for t in range(1, 13)}.isdisjoint({4})

    for r, depth, residue in [(13, 2, 147), (19, 3, 5903)]:
        modulus = r**depth
        t = residue % r
        assert {pow(x, modulus//r, modulus) for x in range(t, modulus, r)} == {residue}
        assert sum(x % modulus for x in (114*residue, 117*residue, 16*residue+1)) >= modulus

    assert 43 % 13 == 4 and 127 % 19 == 13 and pow(127, 18, 361) == 1
    for r, x in [(13, 43), (19, 127)]:
        assert avoids((114*x, 117*x, 16*x+1), r)
    x = Fraction(-16, 121)
    n = 1+247*x
    assert n.numerator % 3 == 0 and n.numerator % 9 != 0 and n.denominator % 3 != 0
    parts = (114*x, 117*x, 16*x+1)
    seen = {}
    while parts not in seen:
        seen[parts] = len(seen)
        digits = tuple(p.numerator*pow(p.denominator, -1, 3) % 3 for p in parts)
        assert sum(digits) < 3
        parts = tuple((p-v)/3 for p, v in zip(parts, digits))
    assert seen[parts] == 4 and len(seen)-seen[parts] == 5


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--exponent", type=int, default=1)
    parser.add_argument("--limit", type=int, default=247)
    args = parser.parse_args()
    if args.exponent < 1 or args.limit < 3:
        parser.error("the exponent must be positive and the cofactor limit at least 3")
    check_247()
    d, shape, count = find_first(args.exponent, args.limit)
    print(f"Scanned {count} candidates at exponent {args.exponent} through first survivor or cofactor {args.limit}.")
    print("First survivor:", d, shape)
