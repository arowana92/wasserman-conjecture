"""Check all-exponent and bounded-exponent theorems at 1+d*2**a for odd d <= 65535."""
from argparse import ArgumentParser
from functools import lru_cache
from math import gcd
from typing import NamedTuple

from carry_blocks import CarryGraph
from search_endpoint import atoms, smallest_factors, splits


SPF = smallest_factors(65535)


@lru_cache(None)
def factors(n):
    return tuple(atoms(n, SPF))


@lru_cache(None)
def scales(n):
    values = [1]
    for prime, power in factors(n):
        if prime == 2:
            continue
        old = values[:]
        k = prime
        while k <= power:
            values.extend(k*v for v in old)
            k *= prime
    return values[1:]


@lru_cache(None)
def assignments(d):
    return [(s, d//s) for s, _, _ in splits(factors(d))]


def companions(u, t):
    # Work with individual coordinates. Distinct coordinates may be the
    # exceptional residue for disjoint sets of primes of t.
    options = {}
    for v in range(1, u-1):
        if all(v % w in (0, u % w) for _, w in factors(t)):
            options[v] = sum(1 << j for j, (_, w) in enumerate(factors(t)) if v % w)
    for i in options:
        if 3*i > u:
            break
        for j in options:
            if j < i:
                continue
            k = u-i-j
            if k < j:
                break
            if k in options and not (options[i]&options[j] or options[j]&options[k]
                                     or options[k]&options[i]):
                yield i, j, k


@lru_cache(None)
def has_companion(u, t):
    return t == 1 or next(companions(u, t), None) is not None


@lru_cache(None)
def logarithms(w):
    table = {}
    r = 1
    while r not in table:
        table[r] = len(table)
        r = 2*r % w
    return table


def progression(d, x, y, u, s):
    targets = [(w, -pow(d, -1, w) % w) for _, w in factors(u)]
    for p, w in factors(s):
        if x % w and y % w:
            return None
        if (x+y) % w:
            if (x+y) % p == 0:
                return None
            targets.append((w, pow(x+y, -1, w)))
    a, h = 0, 1
    for w, r in targets:
        table = logarithms(w)
        if r not in table:
            return None
        b, m = table[r], len(table)
        g = gcd(m, h)
        if (b-a) % g:
            return None
        k = m//g
        a += h*((b-a)//g*pow(h//g, -1, k) % k)
        h *= k
    return a, h


class Shape(NamedTuple):
    d: int
    x: int
    y: int
    z: int
    u: int
    t: int
    a: int
    h: int

    @property
    def slopes(self):
        return self.x, self.y, self.z

    @property
    def forced(self):
        return sorted(factors(self.u) + factors(self.d//self.t))


def templates(limit):
    mask = (1 << limit.bit_length())-1
    for total in range(4, limit):
        if not scales(total):
            continue
        free = mask ^ total
        x = (total-1)&total
        while 2*x > total:
            y = total-x
            admissible = scales(gcd(x, total))
            if admissible:
                upper = max(admissible)**3
                xy = x*y
                z = free
                while z:
                    d = total+z
                    if d <= limit and d % 2:
                        needed = d//gcd(d, xy)
                        cost = 27*needed*needed
                        if cost <= upper and len(factors(d)) >= 2:
                            for u in admissible:
                                if gcd(d, u) != 1 or cost > u**3:
                                    continue
                                for s, t in assignments(d):
                                    if xy % s or 27*t*t > u**3:
                                        continue
                                    if any(2*w+1 > u for _, w in factors(t)):
                                        continue
                                    exponents = progression(d, x, y, u, s)
                                    if exponents and has_companion(u, t):
                                        yield Shape(d, y, x, z, u, t, *exponents)
                    z = (z-1)&free
            x = (x-1)&total


def orbit(a, h, modulus):
    first = pow(2, a, modulus)
    step = pow(2, h, modulus)
    value = first
    while True:
        yield value
        value = value*step % modulus
        if value == first:
            break


def passes_layers(slopes, shift, value, p, modulus):
    # Legendre's factorial floors, evaluated at each power separately.
    parts = [k*value for k in slopes]
    parts[-1] += shift
    total = sum(parts)
    while modulus > 1:
        if total//modulus != sum(v//modulus for v in parts):
            return False
        modulus //= p
    return True


def local_orbit(shape, p, w, companion=None):
    modulus = p
    while modulus*p <= 65536 or (shape.u % p == 0 and modulus <= w):
        modulus *= p
    inverse_u = pow(shape.u, -1, modulus) if companion else None
    for q in orbit(shape.a, shape.h, modulus):
        if shape.u % p == 0 and (shape.d*q+1) % (w*p) == 0:
            continue
        if companion is None:
            if passes_layers(shape.slopes, 1, q, p, modulus):
                return True
        else:
            v = (shape.d*q+1)*inverse_u % modulus
            if passes_layers(companion, 0, v, p, modulus):
                return True
    return False


def passes_orbits(shape):
    if not all(local_orbit(shape, p, w) for p, w in shape.forced):
        return False
    if shape.t == 1:
        return True
    return any(all(local_orbit(shape, p, w, c) for p, w in factors(shape.t))
               for c in companions(shape.u, shape.t))


def live_edges(graph):
    parents = [[] for _ in graph.states]
    live = [i for i, row in enumerate(graph.states) if graph.no_carry(row)]
    seen = set(live)
    for i, edges in enumerate(graph.edges):
        for _, j in edges:
            parents[j].append(i)
    for j in live:
        for i in parents[j]:
            if i not in seen:
                seen.add(i)
                live.append(i)
    return {i: {v: j for v, j in graph.edges[i] if j in seen} for i in seen}


def growth(edges, first_digits, p):
    start = {edges.get(0, {}).get(r) for r in first_digits} - {None}
    states = list(start)
    for i in states:
        for j in edges[i].values():
            if j not in start:
                start.add(j)
                states.append(j)
    if not states:
        return None
    index = {v: i for i, v in enumerate(states)}
    children = [tuple(index[j] for j in edges[i].values()) for i in states]
    counts = [1]*len(states)
    for _ in range(32):
        counts = [sum(counts[j] for j in row) for row in children]
    first = counts
    for _ in range(32):
        counts = [sum(counts[j] for j in row) for row in children]
    # The ratio comparison is exact, including the choice of the maximum.
    num, den = 0, 1
    for a, b in zip(counts, first):
        if a*den > b*num:
            num, den = a, b
    return next(k for k in range(101) if num**100 <= den**100 * p**(32*k))


def check_range(shape, p, edges, limit):
    # Refine the exponent progression so it fixes one more digit of 2^a.
    modulus, stride = 1, shape.h
    levels = []
    while stride <= limit or not levels:
        previous = modulus
        modulus *= p
        unit = pow(2, stride, modulus)
        branches, value = 1, unit
        while value != 1:
            branches += 1
            value = value*unit % modulus
        levels.append((previous, stride, branches))
        stride *= branches
    levels = [(m, s, k, pow(2, s, modulus)) for m, s, k in levels]
    # A short independent exponent orbit can discard most leaves cheaply.
    tables = []
    for r, _ in shape.forced:
        if r == p:
            continue
        m = r
        while m*r <= 65536:
            m *= r
        table = bytes(passes_layers(shape.slopes, 1, v, r, m)
                      for v in orbit(shape.a, shape.h, m))
        tables.append(table)
    precisions = [(r, r**80) for r, _ in shape.forced]
    leaves = 0

    def visit(i, state, e, q):
        nonlocal leaves
        if i == len(levels):
            if e <= 0:
                return
            leaves += 1
            n = (e-shape.a)//shape.h
            if all(table[n % len(table)] for table in tables):
                if all(passes_layers(shape.slopes, 1, pow(2, e, m), r, m)
                       for r, m in precisions):
                    raise ArithmeticError(("unresolved exponent", shape, e))
            return
        previous, step, branches, multiplier = levels[i]
        for _ in range(branches):
            if e > limit:
                break
            next_state = edges[state].get(q//previous % p)
            if next_state is not None:
                visit(i+1, next_state, e, q)
            e += step
            q = q*multiplier % modulus

    visit(0, 0, shape.a, pow(2, shape.a, modulus))
    return leaves


def read_templates(path):
    with open(path) as rows:
        return sorted(Shape(*map(int, row.split())) for row in rows)


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--limit", type=int, default=65535)
    parser.add_argument("--exponents", type=int, default=10**24)
    parser.add_argument("--compare", help="compare with all rows from the C++ enumeration")
    parser.add_argument("--from-templates", help="use a previously generated C++ list")
    args = parser.parse_args()
    if not 1 <= args.limit <= 65535 or args.exponents < 1:
        parser.error("require a cofactor limit in [1,65535] and a positive exponent bound")
    rows = read_templates(args.from_templates) if args.from_templates else sorted(templates(args.limit))
    if args.compare and rows != read_templates(args.compare):
        raise ArithmeticError("template lists differ")
    if any(row.d > args.limit for row in rows):
        raise ValueError("template cofactor exceeds the supplied limit")
    possible = [row for row in rows if passes_orbits(row)]
    print(f"{len(rows)} templates, {len(possible)} pass the local power orbits.", flush=True)
    remaining = []
    leaves = 0
    for row in possible:
        p = row.forced[0][0]
        edges = live_edges(CarryGraph(p, row.slopes, (0, 0, 1)))
        rate = growth(edges, orbit(row.a, row.h, p), p)
        if rate is None:
            continue
        remaining.append(row)
        count = check_range(row, p, edges, args.exponents)
        leaves += count
        print(*row, "base", p, "percent", rate, "leaves", count, flush=True)
    print(f"{len(remaining)} progressions at {len({r.d for r in remaining})} odd cofactors.")
    print(f"{leaves} terminal exponents tested through {args.exponents}, zero survive.")


if __name__ == "__main__":
    main()
