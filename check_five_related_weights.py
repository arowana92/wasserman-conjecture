"""Synchronized carry graphs for the weights 5t+1, 7t+1, 8t+1, 9t+1, 16t+1."""
from argparse import ArgumentParser
from collections import deque
from math import lcm

from carry_blocks import CarryGraph, intersection


class JointGraph:
    """Read one scalar in several carry graphs in the same base."""
    def __init__(self, *graphs):
        self.graphs = graphs
        self.base = graphs[0].base
        assert all(g.base == self.base for g in graphs)
        tables = [[dict(edges) for edges in graph.edges] for graph in graphs]
        initial = (0,)*len(graphs)
        self.states = [initial]
        ids = {initial: 0}
        self.edges = []
        self.accept = 0
        for i, state in enumerate(self.states):
            if all(g.accept >> s & 1 for g, s in zip(graphs, state)):
                self.accept |= 1 << i
            edges = []
            for digit in range(self.base):
                if not all(digit in table[s] for table, s in zip(tables, state)):
                    continue
                nxt = tuple(table[s][digit] for table, s in zip(tables, state))
                if nxt not in ids:
                    ids[nxt] = len(self.states)
                    self.states.append(nxt)
                edges.append((digit, ids[nxt]))
            self.edges.append(edges)

    def accepts_integer(self, t):
        return all(g.accepts_integer(t) for g in self.graphs)

    def suffix_residues(self, nbits, low, modulus):
        """Return high-part residues admitting some accepted word after the fixed bits."""
        state = 0
        for i in range(nbits):
            state = dict(self.edges[state]).get((low >> i) & 1)
            if state is None:
                return set()
        predecessors = [[] for _ in self.states]
        for parent, edges in enumerate(self.edges):
            for bit, child in edges:
                predecessors[child].append((bit, parent))
        reached = {(i, 0) for i in range(len(self.states)) if self.accept >> i & 1}
        todo = deque(reached)
        while todo:
            child, high = todo.popleft()
            for bit, parent in predecessors[child]:
                pair = parent, (2*high+bit) % modulus
                if pair not in reached:
                    reached.add(pair)
                    todo.append(pair)
        return {high for i, high in reached if i == state}


def cut(p, w, a, step=1, initial=0):
    return CarryGraph(p, (a*step, (w-a)*step), (a*initial, (w-a)*initial+1))


def odd_pair(w, a, u, v, b):
    # For b=7 at v=9 write t=7z+3, so (9t+1)/7=9z+4.
    companion = (v, b) == (9, 7)
    step, initial = (7, 3) if companion else (1, 0)
    other = (CarryGraph(2, (9, 18, 36), (4, 8, 16))
             if companion else cut(2, v, b))
    graph = JointGraph(cut(2, w, a, step, initial), other)
    low_t = pow(a, -1, w)
    if v == 7 and low_t % 8 != 3:
        return False
    if v == 9 and ((b == 4 and low_t % 8 != 3) or (b == 8 and low_t % 4 != 3)):
        return False
    modulus = lcm(u, v, 9 if v == 7 else 1)
    low = (low_t-initial)*pow(step, -1, w) % w
    residues = graph.suffix_residues(w.bit_length()-1, low, modulus)
    for high in residues:
        t = step*(w*high+low)+initial
        if (w*t+1) % u or (b*t-1) % v:
            continue
        if v == 7 and t % 9 != 2:
            continue
        return True
    return False


def compare():
    survivors = []
    for w in (8, 16):
        for a in range(3, w, 2):
            for u in range(3, a+1, 2):
                if a % u or (w == 8 and a == 7):
                    continue
                for v, b in ((7, 6), (9, 4), (9, 7), (9, 8)):
                    if odd_pair(w, a, u, v, b):
                        survivors.append((w, a, u, v, b))
    assert survivors == [(16, 3, 3, 7, 6), (16, 11, 11, 9, 4)]
    for a in (3, 5):
        b = a+8
        low = pow(b, -1, 16)
        graph = JointGraph(cut(2, 8, a), cut(2, 16, b))
        # Same a would force a to divide both weights. For the lifted a+8
        # even the necessary divisibility of the first weight by a suffices.
        assert not any((8*(low+16*h)+1) % a == 0
                       for h in graph.suffix_residues(4, low, a))
    print("The only surviving simultaneous assignments are", survivors)


def binary_formula():
    graph = JointGraph(cut(2, 16, 3), cut(2, 7, 6))
    # A component state stores the multiplication carries for each cut.
    def label(i):
        return tuple(v for g, s in zip(graph.graphs, graph.states[i]) for v in g.states[s])
    start = 0
    for bit in (1, 1, 0, 1):
        start = dict(graph.edges[start])[bit]
    names = ((2,9,4,0), (2,11,5,0), (1,5,2,0), (1,4,2,0),
             (2,8,4,0), (0,2,1,0), (0,1,0,0), (0,0,0,0), (1,7,3,0))
    table = (((0,3),(1,1)), ((0,2),), ((1,0),), ((0,5),(1,4)), ((0,3),),
             ((0,6),), ((0,7),(1,8)), ((0,7),), ())
    assert label(start) == names[0]
    for i, row in enumerate(names):
        idx = next(j for j in range(len(graph.states)) if label(j) == row)
        assert [(b, names.index(label(j))) for b,j in graph.edges[idx]] == list(table[i])
        assert bool(graph.accept >> idx & 1) == (i in (0,3,4,5,6,7))
    for m in range(18):
        for k in range(18):
            t = 11+16*5*((8**m-1)//7)+32*8**m*((4**k-1)//3)
            assert graph.accepts_integer(t) and t % 16 == 11
            assert 42*(t+1) == 448*8**m*4**k+32*8**m+24
            assert bin(42*(t+1)).count('1') <= 6
            L, B = 32*8**m, 4**k
            assert 42*t == L*(14*B+1)-18
            assert (t+1) % L + (13*t+1) % L == L-4
            assert ((t+1) & (13*t+1)) == L*((B-1)//3)
    print("The complete binary transition table and the two-exponent formula agree")


def factor_eight_dimensions():
    for a, odd_bound in ((3, 644), (5, 710)):
        for p, numerator in ((2, 488), (a, odd_bound)):
            paths = cut(p, 8, a).block_bound(4096)
            assert paths**1000 < p**(4096*numerator)
            print(f"Factor 8 with a={a}, base {p}: {paths.bit_length()} path-count bits")


def search(digits):
    for v, a, b, modulus, residue in ((7, 3, 6, 1008, 0), (9, 11, 4, 1584, 1267)):
        if v == 7:
            residue = next(t for t in range(modulus) if t % 16 == 11 and t % 72 == 11
                           and 6*t % 7 == 1)
        graphs = [JointGraph(cut(2, 16, a), cut(2, v, b))]
        if v == 7:
            graphs.extend((JointGraph(cut(3, 16, a), cut(3, v, b)), cut(7, v, b)))
        else:
            graphs.extend((cut(3, v, b), cut(11, 16, a)))
        counts = [CarryGraph.block_bound(g, 256) for g in graphs]
        if v == 9:
            for g, count in zip(graphs, counts):
                print("Base", g.base, "256-digit path-count bit length", count.bit_length())
            # The verifier checks upper dimensions for the three block alphabets.
            for g, count, j in zip(graphs, counts, (31, 56, 77)):
                assert count**100 <= g.base**(256*j)
        small, _, _ = intersection(graphs, 80000, modulus, residue)
        expected = [t for t in range(residue, 80001, modulus)
                    if t >= 2 and all(g.accepts_integer(t) for g in graphs)]
        assert small == expected
        ans, nodes, tested = intersection(graphs, 10**digits, modulus, residue)
        print(f"Pair 16 and {v} through 10^{digits}: {nodes:,} nodes, {tested:,} "
              f"direct checks, {len(ans)} survivors", flush=True)
        assert not ans


if __name__ == "__main__":
    cli = ArgumentParser(description=__doc__)
    cli.add_argument("--digits", type=int, default=1000)
    args = cli.parse_args()
    compare()
    binary_formula()
    factor_eight_dimensions()
    search(args.digits)
