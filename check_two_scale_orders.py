"""Certify two-base order exponents at all primitive denominators through 1000."""

from functools import lru_cache
from math import ceil, log

from sympy import Matrix, symbols

from carry_blocks import CarryGraph
from check_logarithmic_intervals import forced_configurations, graph_for


@lru_cache(maxsize=None)
def dimension(base, slopes):
    graph = graph_for(base, slopes)
    predecessors = [[] for _ in graph.states]
    for source, edges in enumerate(graph.edges):
        for _, target in edges:
            predecessors[target].append(source)

    retained = {i for i in range(len(graph.states)) if graph.accept >> i & 1}
    pending = list(retained)
    while pending:
        for source in predecessors[pending.pop()]:
            if source not in retained:
                retained.add(source)
                pending.append(source)

    successors = [
        [target for _, target in edges if target in retained] if i in retained else []
        for i, edges in enumerate(graph.edges)
    ]
    counts = [int(i in retained) for i in range(len(successors))]
    for _ in range(300):
        counts = [sum(counts[k] for k in targets) for targets in successors]
    next_counts = [sum(counts[k] for k in targets) for targets in successors]
    assert all(x or not y for x, y in zip(counts, next_counts))
    if not any(counts):
        return 0

    # Round each ratio upward. The floating-point proposal is checked exactly.
    unit = 10**6
    factor = max((unit*y + x - 1)//x for x, y in zip(counts, next_counts) if x)
    exponent = max(0, ceil(10000*log(factor/unit, base)))
    assert factor**10000 <= unit**10000 * base**exponent
    return exponent


def check():
    results = set()
    cube_shapes = {(119, 203, 664), (203, 714, 69)}
    for denominator, slopes, bases in forced_configurations():
        (s, p), (t, r) = sorted((dimension(b, slopes), b) for b in bases)[:2]
        results.add((denominator, slopes, p, r, s, t))

    worst = max(s+t for _, _, _, _, s, t in results)
    cubes = max(s+t for d, slopes, _, _, s, t in results if d == 986 and slopes in cube_shapes)
    assert len(results) == 16 and (worst, cubes) == (7776, 7203)
    for d, slopes, p, r, s, t in sorted(results):
        print(f'{d} {slopes}: bases {p}, {r}; exponents {s}/10000, {t}/10000.')
    print(f'Largest sum {worst}/10000, at cubes {cubes}/10000.')

    graph = CarryGraph(2, (6, 1), (0, 1))
    matrix = Matrix([
        [sum(k == j for _, k in edges) for j in range(len(graph.states))]
        for edges in graph.edges
    ])
    x = symbols('x')
    expected = x**3 * (x-1) * (x**3-x-1) * (x**3-x**2-1)
    assert matrix.charpoly(x).as_expr().expand() == expected.expand()
    print('Binary multiplier-six characteristic polynomial verified.')


if __name__ == '__main__':
    check()
