"""Check the binary parity criterion and classify the small inverse residues."""

from collections import defaultdict
from math import comb

from carry_blocks import CarryGraph


def factorial_twos(n):
    value = 0
    while n:
        n //= 2
        value += n
    return value


def live_states(graph):
    parents = defaultdict(list)
    for source, edges in enumerate(graph.edges):
        for _, dest in edges:
            parents[dest].append(source)
    # Zero padding eventually reaches state zero precisely for valid endings.
    todo = [0]
    live = {0}
    while todo:
        for state in parents[todo.pop()]:
            if state not in live:
                live.add(state)
                todo.append(state)
    return live


def check_differences(h, b):
    a = 2 ** (h - 2)
    f = (a, 2*a, a - 2**b + 1)
    forbidden = set(range(1, h - b)) | set(range(b, h))
    # Reachable from the start, and with a suffix leading back to zero.
    g = CarryGraph(2, f, (0, 0, 0))
    live = live_states(g)
    window = {t for s in live for digit, t in g.edges[s] if digit == 1 and t in live}
    for distance in range(1, h + 1):
        another_one = any(digit == 1 and t in live for s in window for digit, t in g.edges[s])
        assert another_one == (distance not in forbidden), (h, b, distance)
        window = {t for s in window for _, t in g.edges[s] if t in live}

    for z in range(1, 1025):
        expected = all((z & (z >> gap)) == 0 for gap in forbidden)
        parts = [v*z for v in f]
        v2 = factorial_twos(sum(parts)) - sum(map(factorial_twos, parts))
        assert (v2 == 0) == expected, (h, b, z)
        if sum(parts) <= 500:
            exact = comb(sum(parts), parts[0]) * comb(sum(parts[1:]), parts[1])
            assert (exact % 2 == 1) == expected
    return len(live)


def completion(a, d):
    # Low log2(a) bits of the input are forced. Only d contributes below a.
    prefix = pow(d, -1, a)
    slopes = (a, 2*a, d)
    start = tuple(f*prefix//a for f in slopes)
    graph = CarryGraph(2, slopes, start)
    # Store actual input digits on the first path to each state.
    paths = {0: (0, 0)}
    queue = [0]
    for state in queue:
        suffix, length = paths[state]
        if graph.accept >> state & 1:
            z = prefix + a*suffix
            x, y, t = (f*z for f in slopes)
            assert d*z % a == 1 and (x & y | y & t | x & t) == 0
            return z
        for bit, dest in graph.edges[state]:
            if dest not in paths:
                paths[dest] = (suffix + (bit << length), length + 1)
                queue.append(dest)
    return None


def reverse_completion(a, d):
    # Enumerate every feasible carry state from real prefix fractions, including
    # the total as a fourth coordinate. Build edges backwards by division.
    slopes = (a, 2*a, d, 3*a + d)
    cells = {tuple(f*n//den for f in slopes) for den in slopes for n in range(den)}
    clear = {s for s in cells if sum(s[:3]) == s[3]}
    zero = (0, 0, 0, 0)
    todo = [zero]
    live = {zero}
    units = ((0, 0, 0, 0), (1, 0, 0, 1), (0, 1, 0, 1), (0, 0, 1, 1))
    while todo:
        dest = todo.pop()
        for bit in (0, 1):
            for unit in units:
                src = tuple(2*y + u - bit*f for y, u, f in zip(dest, unit, slopes))
                if src in clear and src not in live:
                    live.add(src)
                    todo.append(src)
    prefix = pow(d, -1, a)
    return tuple(f*prefix//a for f in slopes) in live


def main():
    pairs = states = 0
    for h in range(3, 17):
        for b in range(1, h - 1):
            states += check_differences(h, b)
            pairs += 1
    print(f"{pairs} parameter pairs, {states} live states, {pairs*1024} valuation checks")
    for m in range(2, 11):
        a = 2**m
        found = []
        for d in range(1, a, 2):
            witness = completion(a, d)
            prefix = pow(d, -1, a)
            # Adjacent ones already make the two scaled power-of-two parts overlap.
            possible = (prefix & (prefix >> 1)) == 0 and reverse_completion(a, d)
            assert possible == (witness is not None), (m, d)
            if possible:
                found.append(d)
        print(f"m={m}, admissible d: {found}")


if __name__ == "__main__":
    main()
