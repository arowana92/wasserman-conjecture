"""Finite carry tests for arbitrary arithmetic progressions of binary exponents."""
from collections import Counter
from itertools import combinations
from math import gcd, lcm
from pathlib import Path
from typing import NamedTuple

from carry_blocks import CarryGraph
from check_binary_cofactor_rays import Shape, live_edges, orbit, passes_layers


class Ray(NamedTuple):
    row: Shape
    prime: int
    terms: int


def published_rays():
    # These rows are the output of the exhaustive bounded-cofactor theorem.
    table = Path(__file__).with_name("binary_cofactor_rays.md").read_text()
    four = {33661, 56279, 60829, 61423}
    rays = []
    for line in table.splitlines():
        cells = [v.strip() for v in line.split("|")[1:-1]]
        if len(cells) == 8 and cells[0].isdigit():
            d = int(cells[0])
            slopes = tuple(map(int, cells[1].split(",")))
            row = Shape(d, *slopes, *map(int, cells[2:6]))
            rays.append(Ray(row, int(cells[6]), 4 if d in four else 3))
    if len(rays) != 11:
        raise ValueError("expected eleven published cofactor progressions")
    return rays


def reachable(edges, starts):
    seen = set(starts)
    work = list(seen)
    for i in work:
        for j in edges[i].values():
            if j not in seen:
                seen.add(j)
                work.append(j)
    return seen


def linear_children(edges, p, states, carries, initial):
    # If the first two words are x and y, the ith word is x+i(y-x).
    for a, first in edges[states[0]].items():
        for b, second in edges[states[1]].items():
            if initial and a == b:
                continue
            children, following = [first, second], []
            for i in range(2, len(states)):
                carry, digit = divmod(a+i*(b-a)+carries[i-2], p)
                child = edges[states[i]].get(digit)
                if child is None:
                    break
                following.append(carry)
                children.append(child)
            else:
                yield tuple(children), tuple(following)


def linear_depth(edges, p, states, terms):
    current = {((i,)*terms, (0,)*(terms-2)) for i in states}
    for k in range(1, 33):
        current = {child for previous, carries in current
                   for child in linear_children(edges, p, previous, carries, k == 1)}
        if not current:
            return k
    raise ArithmeticError("linear progression remained possible at the tested depth")


def geometric_depth(edges, p, terms, t, prefix, state):
    # At small t the unit ratio is fixed by the first two complete residues.
    modulus = p**t
    nodes = [(prefix, prefix, (state,)*terms)]
    for depth in range(1, 33):
        following = []
        next_modulus = modulus*p
        for a, b, states in nodes:
            for da, first in edges[states[0]].items():
                x = a+modulus*da
                inverse = pow(x, -1, next_modulus)
                for db, second in edges[states[1]].items():
                    if depth == 1 and da == db:
                        continue
                    y = b+modulus*db
                    value, ratio = y, y*inverse % next_modulus
                    dest = [first, second]
                    for i in states[2:]:
                        value = value*ratio % next_modulus
                        child = edges[i].get(value//modulus)
                        if child is None:
                            break
                        dest.append(child)
                    else:
                        following.append((x, y, tuple(dest)))
        if not following:
            return depth
        nodes, modulus = following, next_modulus
    raise ArithmeticError("geometric progression remained possible at the tested depth")


def lift_progression(row, p0, edges):
    a, h = row.a, row.h
    for p, w in row.forced:
        modulus = p
        while modulus*p <= 65536:
            modulus *= p
        if row.u % p == 0:
            modulus = max(modulus, w*p)
        allowed = []
        for i, q in enumerate(orbit(a, h, modulus)):
            if row.u % p == 0 and (row.d*q+1) % (w*p) == 0:
                continue
            if p == p0:
                state, value, m = 0, q, modulus
                while m > 1 and state is not None:
                    state = edges[state].get(value % p)
                    value //= p
                    m //= p
                if state is None:
                    continue
            elif not passes_layers(row.slopes, 1, q, p, modulus):
                continue
            allowed.append(i)
        if not allowed:
            raise ArithmeticError(("empty local orbit", row, p))
        g = gcd(i+1, *(j-allowed[0] for j in allowed))
        a += h*(allowed[0] % g)
        h *= g
    return row._replace(a=a, h=h)


class Packet:
    def __init__(self, edges, states, p):
        self.p = p
        index = {s: i for i, s in enumerate(states)}
        self.children = [[index[j] for j in edges[s].values()] for s in states]
        self.profile = {0: [0]*len(states), 1: [1]*len(states)}

    def counts(self, k):
        if k not in self.profile:
            small, extra = divmod(k, self.p)
            lo = self.counts(small)
            hi = self.counts(small+1) if extra else lo
            # The p first digits occur in groups of size small or small+1.
            # An extra index increases any subtree optimum by at most one.
            counts = []
            for children in self.children:
                result = sum(lo[j] for j in children)
                increments = sum(hi[j]-lo[j] for j in children)
                counts.append(result + min(extra, increments))
            self.profile[k] = counts
        return self.profile[k]


def greatest_step(rows):
    step = lcm(*(r.h for r in rows))
    for r, s in combinations(rows, 2):
        common = gcd(r.h, s.h)
        difference = r.a-s.a
        while True:
            present = gcd(step, common)
            loss = present//gcd(present, difference)
            if loss == 1:
                break
            step //= loss
    return step


def joint_counts(rows, packets, lengths):
    totals = Counter()
    cache = {}
    for size in range(1, len(rows)+1):
        for indices in combinations(range(len(rows)), size):
            step = greatest_step([rows[i] for i in indices])
            periods = [rows[i].h//gcd(rows[i].h, step) for i in indices]
            for k in lengths:
                upper = 0
                for i, period in zip(indices, periods):
                    n = 1+(k-1)//period
                    if (i, n) not in cache:
                        cache[i, n] = max(packets[i].counts(n))
                    upper += cache[i, n]
                totals[k] = max(totals[k], upper)
    return totals


def main():
    rows, packets = [], []
    for ray in published_rays():
        r, p, terms = ray
        edges = live_edges(CarryGraph(p, r.slopes, (0, 0, 1)))
        first, = orbit(r.a, r.h, p)
        starts = [(first, edges[0][first])]
        states = sorted(reachable(edges, [starts[0][1]]))
        k = linear_depth(edges, p, states, terms)
        prefixes, shallow, checked = starts, 0, 0
        for t in range(1, k):
            shallow = max(shallow, *(geometric_depth(edges, p, terms, t, x, s)
                                     for x, s in prefixes))
            checked += len(prefixes)
            prefixes = [(x+digit*p**t, dest) for x, s in prefixes
                        for digit, dest in edges[s].items()]
        tight = lift_progression(r, p, edges)
        rows.append(tight)
        packets.append(Packet(edges, states, p))
        print(r.d, r.a, "terms", terms, "states", len(states), "depth", k,
              "shallow", checked, shallow, "exponent", tight.a, tight.h, flush=True)
    assert (2*(801-717)) % 1743 != 0
    lengths = (*range(2, 51), 100, 1000, 10**6, 10**9)
    bound = joint_counts(rows, packets, lengths)
    first = next(k for k in lengths if bound[k] < k)
    for k in (first, 50, 100, 1000, 10**6, 10**9):
        print(k, "indices,", bound[k], "possibly exceptional at any cofactor")


if __name__ == "__main__":
    main()
