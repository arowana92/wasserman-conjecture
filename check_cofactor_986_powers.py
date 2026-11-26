"""Compare powers in the two possible rows at 986 using finite carry graphs."""
from argparse import ArgumentParser
from itertools import product
from math import gcd

from carry_blocks import CarryGraph

ROWS = (((119, 203, 664), 4215), ((203, 714, 69), 3501))
INVERSES = (0, 1, 4, 5, 2, 3, 6)


def productive_edges(base, slopes):
    graph = CarryGraph(base, slopes, (0, 0, 1))
    live = {i for i in range(len(graph.states)) if graph.accept >> i & 1}
    while True:
        expanded = live | {i for i, edges in enumerate(graph.edges)
                           if any(t in live for _, t in edges)}
        if expanded == live:
            break
        live = expanded
    return {s: {d: t for d, t in graph.edges[s] if t in live} for s in live}


def low_prefixes(edges, p, first, depth):
    prefixes = [(first, edges[0][first])]
    for j in range(1, depth):
        prefixes = [(x+d*p**j, t) for x, s in prefixes for d, t in edges[s].items()]
    return prefixes


def residue_checks():
    for (slopes, R), expected in zip(ROWS, ({22, 29}, {22})):
        for p, values in ((7, expected), (17, {84})):
            edges = productive_edges(p, slopes)
            assert {x for x, _ in low_prefixes(edges,p,R%p,2)} == values
    assert pow(10,19,29) == 21
    assert {pow(10,i,29) for i in range(28)} == set(range(1,29))
    assert pow(19,-1,28) == 3 and not {19%7, 3%7} & {1,6}
    single, = low_prefixes(productive_edges(7,ROWS[1][0]),7,1,5)
    assert single[0] == 5608
    print("Second row, unique prefix modulo 7^5:", single[0])


def packet(slopes, limit):
    edges = productive_edges(7, slopes)
    root = edges[0][1]
    nodes = depth = 0
    for residue in range(1, 476):
        if gcd(residue,476) != 1:
            continue
        candidates = [(m,root) for m in range(residue,limit+1,476)]
        if len(candidates) < 2:
            continue
        stack = [(1,7,1,candidates)]
        while stack:
            k, modulus, x, candidates = stack.pop()
            assert k < 64, ("two indices unresolved",x,candidates)
            nodes += 1
            depth = max(depth,k)
            following = [[] for _ in range(7)]
            next_modulus = 7*modulus
            for m, state in candidates:
                constant = pow(x,m,next_modulus)//modulus
                inverse = INVERSES[m%7]
                for output_digit, target in edges[state].items():
                    digit = (output_digit-constant)*inverse % 7
                    following[digit].append((m,target))
            for digit, group in enumerate(following):
                if len(group) >= 2:
                    stack.append((k+1,next_modulus,x+digit*modulus,group))
    print(slopes, "packet limit", limit, "nodes", nodes, "maximum depth", depth, flush=True)


def geometric_prefixes(edges, r, state, h, terms, depth_limit):
    """Increment W0,W1 and reconstruct the geometric progression modulo each layer."""
    modulus = 7**h
    nodes = [(r,r,(state,)*terms)]
    for k in range(1, depth_limit+1):
        following = []
        next_modulus = 7*modulus
        for a,b,states in nodes:
            for da,sa in edges[states[0]].items():
                aa = a+modulus*da
                inverse = pow(aa,-1,next_modulus)
                for db,sb in edges[states[1]].items():
                    if k == 1 and da == db:
                        continue
                    bb = b+modulus*db
                    ratio = bb*inverse % next_modulus
                    value = bb
                    new_states = [sa,sb]
                    for s in states[2:]:
                        value = value*ratio % next_modulus
                        digit = value//modulus
                        if digit not in edges[s]:
                            break
                        new_states.append(edges[s][digit])
                    if len(new_states) == terms:
                        following.append((aa,bb,tuple(new_states)))
        if not following:
            return k
        nodes = following
        modulus = next_modulus
    raise AssertionError(("geometric progression unresolved",r,h,terms))


def progression_checks(slopes):
    edges = productive_edges(7,slopes)
    root = edges[0][1]
    reach = {root}
    while True:
        expanded = reach | {t for s in reach for t in edges[s].values()}
        if expanded == reach:
            break
        reach = expanded

    # At h >= 2 the normalized terms, modulo 49, are an ordinary progression.
    for state in reach:
        words = {d+7*e for d,t in edges[state].items() for e in edges[t]}
        for a,b in product(words,repeat=2):
            if (a-b)%7:
                assert (2*b-a)%49 not in words or (3*b-2*a)%49 not in words
    geometric_prefixes(edges,1,root,1,4,2)

    prefixes = [(1,root)]
    peaks = []
    for h in range(1,19):
        peaks.append(max(geometric_prefixes(edges,r,s,h,3,64) for r,s in prefixes))
        prefixes = [(r+7**h*d,t) for r,s in prefixes for d,t in edges[s].items()]
    print(slopes,"three-term lift depths for h=1..18:",peaks,flush=True)


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--limit',type=int,default=10**6)
    args = parser.parse_args()
    residue_checks()
    for slopes, _ in ROWS:
        progression_checks(slopes)
        packet(slopes,args.limit)
