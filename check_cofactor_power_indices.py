"""Classify power indices by local unit powers for cofactors through 1000 and at 2278."""

from collections import Counter, deque
from fractions import Fraction
from math import gcd, lcm

from check_higher_successor_ray import finishing_graph, read_prefix
from check_power_obstruction_cover import DENOMINATORS
from check_prime_power_cofactors import shapes
from check_shifted_directions import possible_scales
from check_shifted_directions_cross import independent_shapes
from search_endpoint import atoms, avoids, smallest_factors


EXPECTED = (2, 21, 51, 595, 2523, 17017, 24871, 40375, 78625, 83657, 272935, 382109,
            531505, 2347835, 8921773, 42632447, 67148081, 166654111, 2573653511)
EXTRA = (3, 8, 49, 67, 289, 2756)


def valuation(n, r):
    if n < 1:
        raise ValueError("valuation needs a positive integer")
    k = 0
    while n % r == 0:
        n //= r
        k += 1
    return k


def records():
    # For 2278 exhaust every divisor that could be the primitive denominator.
    divisors = {m for m in range(3, 2279) if 2278 % m == 0}
    spf = smallest_factors(2278)
    found = []
    for m in sorted(divisors | set(DENOMINATORS)):
        factors = atoms(m, spf)
        rows = {s for s in shapes(m, factors) if gcd(*s[1:4]) == 1}
        if m in divisors and rows != independent_shapes(m):
            raise ArithmeticError(("independent shape enumeration", m))
        for shape in sorted(rows):
            for parity, prog, first in possible_scales(m, factors, shape, spf):
                if m in divisors and (m != 2278 or parity != 1):
                    raise ArithmeticError(("additional 2278 assignment", m, shape, parity))
                if parity:
                    found.append((m, shape, prog, first))
    assert len(found) == 16 and sum(d == 2278 for d, *_ in found) == 2
    return found


def finish(r, slopes, edges, state, low, w):
    """Return a finite input continuing this prefix to a carry-free row."""
    todo = deque([(low, w, state)])
    seen = {state}
    x, y, z = slopes
    while todo:
        v, place, state = todo.popleft()
        if avoids((x*v, y*v, z*v+1), r):
            return v
        for digit, child in edges[state].items():
            if child not in seen:
                seen.add(child)
                todo.append((v+digit*place, place*r, child))
    return None


def compatible_prefix(d, u, slopes, R, M, r, root, w, edges):
    base_part = r**valuation(M, r)
    full_u = r**valuation(u, r)
    precision = max(w, base_part, r*full_u if full_u > 1 else 1)
    for v in range(root, precision, w):
        if (v-R) % base_part or (full_u > 1 and (d*v+1) % (r*full_u) == 0):
            continue
        state = read_prefix(edges, r, v, valuation(precision, r))
        if state is not None:
            witness = finish(r, slopes, edges, state, v, precision)
            if witness is not None:
                return witness
    return None


def local_obstructions(d, shape, prog, r):
    u = shape[0]
    slopes = shape[1:4]
    R, M = prog
    t = R % r
    assert M % r == 0 and t != 0
    order = next(j for j in range(1, r) if pow(t, j, r) == 1)
    obstructions = {s**(valuation(r-1, s)-valuation(order, s)+1)
                    for s, _ in atoms(order, smallest_factors(order))}
    edges = finishing_graph(r, slopes, (0, 0, 1))
    root, w = t, r
    last = None
    while True:
        candidate = compatible_prefix(d, u, slopes, R, M, r, root, w, edges)
        if candidate is None:
            if last is None:
                raise ArithmeticError(("infeasible listed shape", d, slopes, r))
            obstructions.add(w//r)
            return obstructions, last, (root, w)
        if candidate == root == 1:
            return obstructions, 1, None
        last = candidate
        w *= r
        root = pow(root, r, w)


def minimal(values):
    kept = []
    for n in sorted(set(values)):
        if all(n % k for k in kept):
            kept.append(n)
    return tuple(kept)


def common_indices(rows):
    covered = (1,)
    for row in rows:
        covered = minimal(lcm(a, b) for a in covered for b in row)
    return covered


def density(indices):
    coefficients = Counter()
    for n in indices:
        for m, c in list(coefficients.items()):
            coefficients[lcm(m, n)] -= c
        coefficients[n] += 1
        coefficients = Counter({m: c for m, c in coefficients.items() if c})
    return sum(Fraction(c, m) for m, c in coefficients.items())


def solve():
    small, extra = [], []
    for d, shape, prog, first in records():
        # The earlier theorem settles even exponents throughout the first 1000.
        barred = {2} if d <= 1000 else {2756}
        witnesses = {}
        for r in first:
            if r != 2:
                values, witness, _ = local_obstructions(d, shape, prog, r)
                barred.update(values)
                witnesses[r] = witness
        result = minimal(barred)
        print(d, shape[1], result, witnesses, flush=True)
        (small if d <= 1000 else extra).append(result)
    assert common_indices(small) == EXPECTED and common_indices(extra) == EXTRA
    for indices in (EXPECTED, EXTRA):
        print("Minimal common indices", indices, "density", density(indices))


if __name__ == "__main__":
    solve()
