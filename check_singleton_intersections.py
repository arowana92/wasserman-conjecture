"""Check simultaneous multiplicative carries and dimension bounds for the remaining shapes."""

from argparse import ArgumentParser
from collections import Counter
from concurrent.futures import ProcessPoolExecutor

from carry_blocks import CarryGraph, intersection
from check_shifted_prefixes import floor_cover
from check_singleton_scales import ZERO, atoms, surviving_shapes


def percent(graph, length=128):
    # An accepted digit block must belong to a path that can reach acceptance.
    rev = [set() for _ in graph.states]
    for i, edges in enumerate(graph.edges):
        for _,j in edges:
            rev[j].add(i)
    good = {i for i in range(len(rev)) if (graph.accept >> i) & 1}
    stack = list(good)
    while stack:
        new = rev[stack.pop()] - good
        good.update(new)
        stack.extend(new)
    ways = [int(i in good) for i in range(len(rev))]
    for _ in range(length):
        ways = [sum(ways[j] for _,j in edges) if i in good else 0
                for i,edges in enumerate(graph.edges)]
    total = sum(ways)**100
    return next(j for j in range(151) if total <= graph.base**(j*length))


def check(data):
    shape,digits,floor_digits = data
    U,T,k,f,g,y = shape
    bases = [q for q,_ in atoms(T)]
    graphs = [CarryGraph(q,y,ZERO) for q in bases]
    exponents = sorted(percent(graph) for graph in graphs)
    pair_bound = sum(exponents[:2])
    if pair_bound > 100:
        # The first row then avoids 2 in the two shapes requiring this argument.
        assert U%2 == 0 and 2 not in bases
        # The binary affine three-part bound is below 56/100.
        assert 56 + exponents[0] <= 100
    residue = pow(U,-1,T)
    survivors, nodes, _ = intersection(graphs,10**digits,T,residue)
    assert not survivors, (shape,survivors)
    survivors, steps = floor_cover(y,bases,T,residue,10**floor_digits,offsets=ZERO)
    assert not survivors, (shape,survivors)
    return nodes, steps, pair_bound > 100


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--digits', type=int, default=1000)
    parser.add_argument('--floor-digits', type=int, default=50)
    parser.add_argument('--jobs', type=int, default=1)
    args = parser.parse_args()
    assert min(args.digits,args.floor_digits,args.jobs) > 0
    # For lambda^3=lambda^2+1, lambda<147/100<2^(56/100).
    assert 147**3 > 147**2*100+100**3 and 147**100 < 100**100*2**56
    shapes = list(surviving_shapes(5000,Counter()))
    with ProcessPoolExecutor(args.jobs) as pool:
        result = list(pool.map(check,[(s,args.digits,args.floor_digits) for s in shapes]))
    print(f'{len(result)} shapes, {sum(r[0] for r in result)} carry-prefix nodes, '
          f'{sum(r[1] for r in result)} factorial-floor intervals, '
          f'{sum(r[2] for r in result)} dimension comparisons using both rows.')


if __name__ == '__main__':
    main()
