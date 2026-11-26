"""Classify the base-13 carry language of the (104,1524,23) row at cofactor 1651."""

from collections import defaultdict

from carry_blocks import CarryGraph


SLOPES = (104,1524,23)


def check_graph():
    graph = CarryGraph(13,SLOPES,(0,0,1))
    initial = dict(graph.edges[0])[2297%13]
    reachable = {initial}
    pending = [initial]
    while pending:
        for _,j in graph.edges[pending.pop()]:
            if j not in reachable:
                reachable.add(j)
                pending.append(j)
    reverse = defaultdict(set)
    for i in reachable:
        for _,j in graph.edges[i]:
            reverse[j].add(i)

    def ancestors(seeds):
        found, todo = set(seeds), list(seeds)
        while todo:
            for i in reverse[todo.pop()]:
                if i not in found:
                    found.add(i)
                    todo.append(i)
        return found

    live = ancestors({i for i in reachable if (graph.accept >> i) & 1})
    edges = {i: [(d,j) for d,j in graph.edges[i] if j in live] for i in live}
    periodic = graph.states.index((69,1016,15))
    zero = graph.states.index((0,0,0))
    assert {(i,d) for i in live for d,j in edges[i] if j == i} == {(periodic,8),(zero,0)}
    assert edges[zero] == [(0,zero)]
    # Topological deletion checks that these are the only cycles.
    rest = {i:{j for _,j in edges[i] if j != i} for i in live}
    longest = {}
    while rest:
        ends = [i for i,js in rest.items() if js <= longest.keys()]
        assert ends
        for i in ends:
            longest[i] = max((1+longest[j] for j in rest.pop(i)), default=0)
    assert longest[initial] == 72
    after = {periodic}
    pending = [periodic]
    while pending:
        i = pending.pop()
        for d,j in edges[i]:
            if i == j == periodic:
                continue
            assert d == 0
            if j not in after:
                after.add(j)
                pending.append(j)

    toward = ancestors({periodic})

    def prefixes(i, length, value):
        if i == periodic:
            yield length, 3*value-2*13**length
        else:
            for d,j in edges[i]:
                if j in toward:
                    yield from prefixes(j,length+1,value+d*13**length)

    families = sorted(prefixes(initial,1,2297%13))
    assert len(families) == 8
    survivors = []
    modulus = 3*6604
    for length, c in families:
        power = pow(13,length,modulus)
        assert pow(13,63,modulus)*power%modulus == power
        for n in range(length,length+63):
            if (2*pow(13,n,modulus)+c-3*2297)%modulus == 0:
                survivors.append((length,c,n%63))
    assert survivors == [(3,1,54)]
    # A canonical expansion which visits neither nonzero loop nor any
    # repeated padding has length <= 73. The existing dual search covers it.
    assert 13**74 < 10**1000
    print(len(graph.states), 'states;', len(live), 'live after the required units digit;',
          len(families), 'prefix families;', 'survivor', survivors[0])

    for n in range(4):
        Q = (2*13**n+1)//3
        if n == 3:
            assert Q%5 == 0
        else:
            assert sum((a*Q+b)%32 for a,b in zip(SLOPES,(0,0,1))) >= 32 or (
                sum((a*Q+b)%16 for a,b in zip(SLOPES,(0,0,1))) >= 16)
    assert (2*(13**4-1))%32 == 0 and (13**4-1)%5 == 0
    assert sum(a%5 for a in SLOPES) >= 5


if __name__ == '__main__':
    check_graph()
