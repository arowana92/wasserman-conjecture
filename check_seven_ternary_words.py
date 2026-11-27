"""Count arbitrary internal ternary blocks in each residual companion addition."""

from collections import defaultdict
from math import lcm

from check_seven_quadratic_suffixes import RESIDUAL


def no_carry(numbers):
    if min(numbers) < 0:
        return False
    while any(numbers):
        if sum(x%3 for x in numbers) >= 3:
            return False
        numbers = tuple(x//3 for x in numbers)
    return True


def offsets(D, z, w):
    slopes = (2*z,2*w,D-2*z-2*w)
    common = lcm(*slopes)
    cuts = {k for a in slopes for k in range(0,common,common//a)}
    states = sorted({(slopes[0]*k//common,slopes[1]*k//common,slopes[2]*k//common-e)
                     for k in cuts for e in (0,1)})
    return slopes, states


def graph_data(D, z, w):
    slopes,states = offsets(D,z,w)
    ids = {s:i for i,s in enumerate(states)}
    edges = [[] for s in states]
    reverse = [defaultdict(list) for _ in range(3)]
    terminal = set()
    for i,s in enumerate(states):
        if no_carry(s):
            terminal.add(i)
        for d in range(3):
            values = tuple(x+d*a for x,a in zip(s,slopes))
            if sum(x%3 for x in values) < 3:
                dest = ids[tuple(x//3 for x in values)]
                edges[i].append(dest)
                reverse[d][dest].append(i)
    live, queue = set(terminal), list(terminal)
    for i in queue:
        for table in reverse:
            for old in table[i]:
                if old not in live:
                    live.add(old)
                    queue.append(old)
    ids_in_live = {i:j for j,i in enumerate(sorted(live))}
    internal = [[ids_in_live[v] for v in edges[i] if v in live] for i in sorted(live)]
    return reverse,terminal,internal


def leading_graph(reverse,terminal):
    groups,known = [],{}
    def add(items):
        s = frozenset(items)
        if not s:
            return -1
        if s not in known:
            known[s] = len(groups)
            groups.append(s)
        return known[s]

    def before(s,d):
        return add(v for j in s for v in reverse[d][j])

    first = [(d,i) for d in (1,2) if (i := before(terminal,d)) >= 0]
    transitions = []
    for s in groups:
        transitions.append([(d,i) for d in range(3) if (i := before(s,d)) >= 0])
    return first,transitions


def certify():
    expected = {4587:(2498,354,1115844),7383:(3304,654,1107472),9453:(3550,759,1155459)}
    for D,_,_,z,w in RESIDUAL:
        rev,terminal,edges = graph_data(D,z,w)
        first,transitions = leading_graph(rev,terminal)
        counts = [1]*len(edges)
        for _ in range(512):
            counts = [sum(counts[t] for t in row) for row in edges]
        nxt = [sum(counts[t] for t in row) for row in edges]
        assert all(x or not y for x,y in zip(counts,nxt))
        ratio = max((10**6*y+x-1)//x for x,y in zip(counts,nxt) if x)
        assert (len(edges),len(transitions),ratio) == expected[D]
        assert all(10**6*y <= ratio*x for x,y in zip(counts,nxt))
        assert ratio**38 < 3**5 * (10**6)**38
        print(f'D={D}: {len(edges)} coaccessible offsets, {len(transitions)} leading states, '
              f'one-step path ratio at depth 512 at most {ratio}/1000000.')


if __name__ == '__main__':
    certify()
