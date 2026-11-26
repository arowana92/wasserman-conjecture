"""Check all exponents in 1+2^e P Q for two bounded odd prime powers."""

from argparse import ArgumentParser
from collections import Counter, deque
from functools import lru_cache
from itertools import combinations
from math import gcd
from typing import NamedTuple

from carry_blocks import CarryGraph
from check_cofactor_digit_gap import atoms
from check_three_factors import powers


class Shape(NamedTuple):
    P: int
    Q: int
    k: int
    l: int
    u: int
    first: int
    period: int

    @property
    def slopes(self):
        x, y = self.k*self.Q, self.l*self.P
        return x, y, self.P*self.Q-x-y

    @property
    def bases(self):
        return self.u, next(iter(powers(self.P))), next(iter(powers(self.Q)))


@lru_cache(None)
def log_table(modulus):
    logs = {}
    value = 1
    while value not in logs:
        logs[value] = len(logs)
        value = 2*value % modulus
    return logs


def intersect(first, second):
    a, m = first
    b, n = second
    g = gcd(m, n)
    if (b-a) % g:
        return None
    n //= g
    step = (b-a)//g * pow(m//g, -1, n) % n
    return a+m*step, m*n


def multipliers(P, Q):
    logs = log_table(P)
    d = P*Q
    for k in range(1, P):
        x = k*Q
        if (x & d) != x or gcd(x, P) != 1:
            continue
        e = logs.get(pow(x, -1, P))
        if e is not None:
            yield k, x, (e, len(logs))


def generate(bound):
    counts = Counter()
    shapes = []
    odd_powers = sorted((w, p) for p, w in atoms(bound) if p != 2)
    for (P, p), (Q, q) in combinations(odd_powers, 2):
        if p == q:
            continue
        counts['pairs'] += 1
        d = P*Q
        left = list(multipliers(P, Q))
        right = list(multipliers(Q, P)) if left else []
        for k, x, xclass in left:
            for l, y, yclass in right:
                common = gcd(k, l)
                if x & y or x+y >= d or common < 3:
                    continue
                joint = intersect(xclass, yclass)
                if joint is None:
                    continue
                counts['binary rows with a compatible predecessor orbit'] += 1
                for u in powers(common):
                    if u == 2 or d % u == 0:
                        continue
                    logs = log_table(u)
                    e = logs.get(-pow(d, -1, u) % u)
                    if e is not None:
                        full = intersect(joint, (e, len(logs)))
                        if full is not None:
                            shapes.append(Shape(P, Q, k, l, u, *full))
    return counts, shapes


def fixed_modulus(shape, p):
    modulus = p
    assert pow(2, shape.period, p) == 1
    while pow(2, shape.period, modulus*p) == 1:
        modulus *= p
    return modulus


def first_empty_level(shape, p, levels, budget):
    fixed = fixed_modulus(shape, p)
    target = pow(2, shape.first, fixed)
    x, y, z = shape.slopes
    d, common = shape.P*shape.Q, gcd(shape.k, shape.l)
    residues = [0]
    old = 1
    for depth in range(1, levels+1):
        modulus = old*p
        if modulus <= fixed:
            candidates = (target % modulus,)
        elif p*len(residues) <= budget:
            candidates = (v+j*old for v in residues for j in range(p))
        else:
            break
        residues = [t for t in candidates
                    if x*t % modulus + y*t % modulus + (z*t+1) % modulus < modulus
                    and (p != shape.u or (d*t+1) % modulus or common % modulus == 0)]
        if not residues:
            return depth
        old = modulus
    return None


def terminal_possible(shape):
    p = shape.u
    fixed = fixed_modulus(shape, p)
    t = pow(2, shape.first, fixed)
    x, y, z = shape.slopes
    start = (x*t//fixed, y*t//fixed, (z*t+1)//fixed)
    queue = deque([start])
    seen = {start}
    while queue:
        state = queue.popleft()
        row = state
        while any(row) and sum(v % p for v in row) < p:
            row = tuple(v//p for v in row)
        if not any(row):
            return True
        for digit in range(p):
            parts = tuple(v + digit*a for v, a in zip(state, (x, y, z)))
            if sum(v % p for v in parts) >= p:
                continue
            dest = tuple(v//p for v in parts)
            if dest not in seen:
                seen.add(dest)
                queue.append(dest)
    return False


def unlimited(bound):
    counts, shapes = generate(bound)
    counts['scale prime and initial exponent progression'] = len(shapes)
    remaining = []
    for shape in shapes:
        if any(first_empty_level(shape, p, 24, 200000) for p in shape.bases):
            continue
        counts['after local residues'] += 1
        if terminal_possible(shape):
            remaining.append(shape)
    return counts, remaining


def predecessor_table(shape, power):
    base = next(iter(powers(power)))
    modulus = power*base
    x, y, z = shape.slopes
    start = pow(2, shape.first, modulus)
    step = pow(2, shape.period, modulus)
    value = start
    table = []
    while not table or value != start:
        table.append(x*value % modulus + y*value % modulus + (z*value+1) % modulus < modulus)
        value = value*step % modulus
    return table


def ternary_companions(shape):
    # At this pair the odd part of the common multiplier is exactly 47.
    common = gcd(shape.k, shape.l)
    while common % 2 == 0:
        common //= 2
    assert common == shape.u == 47
    p = 3
    w = fixed_modulus(shape, p)
    residue = (shape.P*shape.Q*pow(2, shape.first, w)+1)*pow(47, -1, w) % w
    assert (w, residue) == (81, 51)
    companions = []
    for f in range(1, 47//3+1):
        for g in range(f, (47-f)//2+1):
            row = f, g, 47-f-g
            if any(sum(a*residue % m for a in row) >= m for m in (3, 9, 27, 81)):
                continue
            graph = CarryGraph(p, row, tuple(a*residue//w for a in row))
            if graph.accept:
                companions.append((row, graph))
    assert [row for row, _ in companions] == [(8, 15, 24)]
    return companions[0]


def companion_check(shape, limit, tables):
    row, graph = ternary_companions(shape)
    # Keep precisely the edges from which a finite scalar can terminate.
    live = {i for i in range(len(graph.states)) if graph.accept >> i & 1}
    while True:
        extended = live | {i for i, edges in enumerate(graph.edges)
                           if any(dest in live for _, dest in edges)}
        if extended == live:
            break
        live = extended
    edges = [{digit: dest for digit, dest in choices if dest in live}
             for choices in graph.edges]
    maximum = (limit-shape.first)//shape.period
    stride, precision = 1, 81
    levels = []
    while stride <= maximum:
        precision *= 3
        levels.append((precision, stride))
        stride *= 3
    inverse = pow(47, -1, precision)
    levels = [(w, step, pow(2, shape.period*step, precision)) for w, step in levels]
    leaves = filtered = 0
    x, y, z = shape.slopes
    d = shape.P*shape.Q

    def visit(depth, index, value, state):
        nonlocal leaves, filtered
        if depth == len(levels):
            exponent = shape.first+shape.period*index
            if not 1 <= exponent <= limit:
                return
            leaves += 1
            if any(not table[index % len(table)] for table in tables):
                return
            filtered += 1
            for p in shape.bases:
                t = pow(2, exponent, p**64)
                a, b, c = x*t, y*t, z*t+1
                w = p
                for _ in range(64):
                    if a % w + b % w + c % w >= w:
                        return
                    w *= p
            raise ArithmeticError(('requires more carry tests', shape, exponent, row))
        w, step, multiplier = levels[depth]
        for _ in range(3):
            if index > maximum:
                break
            digit = ((d*value+1)*inverse % w)//(w//3)
            dest = edges[state].get(digit)
            if dest is not None:
                visit(depth+1, index, value, dest)
            index += step
            value = value*multiplier % precision

    visit(0, 0, pow(2, shape.first, precision), 0)
    return leaves, filtered


def finite_check(shape, limit, use_companion=True):
    # When at most two classes survive a whole predecessor digit, intersect
    # those classes before lifting the scale prime. The resulting exponents
    # are disjoint subprogressions of the original one.
    tables = [predecessor_table(shape, power) for power in (shape.P, shape.Q)]
    sparse = next((table for table in tables if len(table) > 1 and sum(table) <= 2), None)
    if sparse is None:
        if use_companion and (shape.P, shape.Q, shape.u) == (379, 797, 47):
            return companion_check(shape, limit, tables)
        return finite_progression(shape, limit, tables)
    counts = [0, 0]
    for j, keep in enumerate(sparse):
        if keep:
            sub = shape._replace(first=shape.first+j*shape.period, period=shape.period*len(sparse))
            a, b = finite_check(sub, limit, use_companion)
            counts[0] += a
            counts[1] += b
    return tuple(counts)


def finite_progression(shape, limit, tables):
    u = shape.u
    x, y, z = shape.slopes
    d = shape.P*shape.Q
    max_index = (limit-shape.first)//shape.period
    # Write e=first+period*index. Each additional scale-prime digit fixes
    # one more base-u digit of index once the initially fixed digits end.
    stride, precision = 1, 1
    levels = []
    while stride <= max_index:
        precision *= u
        branches = 1 if pow(2, shape.period*stride, precision) == 1 else u
        levels.append((precision, stride, branches))
        stride *= branches
    # Carry extra precision while lifting. Each individual exponent can then
    # be screened at eight further layers without modular exponentiation.
    extra_layers = [precision*u**i for i in range(1, 9)]
    tracking = extra_layers[-1]
    digit_steps = []
    for w, step, branches in levels:
        multiplier = pow(2, shape.period*step, tracking)
        if branches == 1:
            digit_steps.append((w, step, (1,), 0))
            continue
        old = w//u
        delta = z*pow(2, shape.first, u)*((multiplier % w - 1)//old) % u
        shifts = tuple(pow(multiplier, j, tracking) for j in range(u))
        digit_steps.append((w, step, shifts, pow(delta, -1, u)))
    leaves = filtered = 0

    def visit(depth, index, value):
        nonlocal leaves, filtered
        if depth == len(digit_steps):
            exponent = shape.first + shape.period*index
            if not 1 <= exponent <= limit:
                return
            leaves += 1
            if any(not table[index % len(table)] for table in tables):
                return
            filtered += 1
            a, b, c = x*value, y*value, z*value+1
            for w in extra_layers:
                if a % w + b % w + c % w >= w:
                    return
            for p in shape.bases:
                residue = pow(2, exponent, p**64)
                a, b, c = x*residue, y*residue, z*residue+1
                w = p
                for _ in range(64):
                    if a % w + b % w + c % w >= w:
                        return
                    w *= p
            raise ArithmeticError(('requires more carry tests', shape, exponent))
        w, step, shifts, inverse_delta = digit_steps[depth]
        subtotal = x*value % w + y*value % w
        if subtotal >= w:
            return
        third = (z*value+1) % w
        if len(shifts) == 1:
            if subtotal+third < w:
                visit(depth+1, index, value)
            return
        # x and y are multiples of u, so their residues are constant over
        # the u lifts. The top digit of z*T+1 takes each value once.
        old = w//u
        top, low = divmod(third, old)
        maximum = (w-1-subtotal-low)//old
        for digit in range(maximum+1):
            branch = (digit-top)*inverse_delta % u
            following = index+step*branch
            if following <= max_index:
                visit(depth+1, following, value*shifts[branch] % tracking)

    visit(0, 0, pow(2, shape.first, tracking))
    return leaves, filtered


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--bound', type=int, default=1000)
    parser.add_argument('--exponent', type=int, default=10**24)
    parser.add_argument('--scale-only', action='store_true',
                        help='also use the first row alone at the pair 379,797')
    args = parser.parse_args()
    if args.bound < 3 or args.exponent < 1:
        parser.error('bounds must be positive, with odd power bound at least 3')
    counts, remaining = unlimited(args.bound)
    print(dict(counts), flush=True)
    print('Unrestricted pairs still possible:', sorted({(s.P, s.Q) for s in remaining}), flush=True)
    totals = [0, 0]
    for s in remaining:
        count = finite_check(s, args.exponent, not args.scale_only)
        totals = [a+b for a, b in zip(totals, count)]
        print(s, 'finite leaves and filtered leaves', count, flush=True)
    print(f'Checked all exponents through {args.exponent}, leaf totals {totals}.')
