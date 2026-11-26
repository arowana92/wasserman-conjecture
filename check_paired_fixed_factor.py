"""Verify the complete accepted binary words for two simultaneous exceptions."""
from argparse import ArgumentParser
from collections import Counter, defaultdict
from math import gcd, isqrt
from pathlib import Path
from subprocess import run


def normalized_rows(p, v, c, e):
    if v == 7:
        return (6*p,c,p-c), (3*p,e,13*p-e)
    return (4*p,c,5*p-c), (11*p,c,5*p-c)


def accepted_inputs(p, v, c, e):
    rows = normalized_rows(p,v,c,e)
    low = ((11 if v == 7 else 3)*pow(p,-1,16)) % 16
    # Starting independently from the four fixed low bits, multiply out directly.
    values = tuple(a*low + (i == 2) for row in rows for i,a in enumerate(row))
    for i in (0,3):
        if any((values[j] & values[k] & 15)
                for j in range(i,i+3) for k in range(j+1,i+3)):
            return []
    initial = tuple(x//16 for x in values)
    states = [initial]
    index = {initial:0}
    edges = []
    parents = defaultdict(list)
    stopping = set()
    slopes = rows[0]+rows[1]
    for i,state in enumerate(states):
        if all(sum(group) == (group[0] | group[1] | group[2])
               for group in (state[:3], state[3:])):
            stopping.add(i)
        edges.append([])
        for bit in (0,1):
            values = tuple(carry+a*bit for carry,a in zip(state,slopes))
            if any(sum(x%2 for x in group)>1 for group in (values[:3],values[3:])):
                continue
            child = tuple(x//2 for x in values)
            if child not in index:
                index[child] = len(states)
                states.append(child)
            j = index[child]
            edges[i].append((bit,j))
            parents[j].append(i)

    # A nonempty canonical suffix ends at its most significant 1. Delete states
    # unable to reach such an edge; this removes arbitrary terminal zero padding.
    productive = {i for i,row in enumerate(edges)
                  if any(bit == 1 and j in stopping for bit,j in row)}
    todo = list(productive)
    while todo:
        for i in parents[todo.pop()]:
            if i not in productive:
                productive.add(i)
                todo.append(i)
    found = [low] if 0 in stopping else []

    def visit(i, value, place, ancestors):
        if i not in productive:
            return
        if i in ancestors:
            raise ValueError(f"Unresolved infinite language at {(p,v,c,e)}")
        for bit,j in edges[i]:
            following = value + bit*place
            if bit and j in stopping:
                found.append(following)
            visit(j, following, 2*place, ancestors | {i})

    visit(0,low,16,set())
    return sorted(found)


def small_powers(limit):
    """Trial-divide each integer, independently of the C++ sieve."""
    for p in range(3,limit+1,2):
        base = next((i for i in range(3, isqrt(p)+1, 2) if p%i == 0), p)
        remaining = p
        while remaining % base == 0:
            remaining //= base
        if remaining == 1:
            yield p


def small_shapes(limit):
    for p in small_powers(limit):
        if gcd(p,42) == 1:
            for c in range(42,p,42):
                for e in (c,c+12*p):
                    if e%48 == 0 and gcd(c,p) == 1:
                        yield p,7,c,e
        if gcd(p,1584) == 1:
            for c in range(1584,5*p,1584):
                if gcd(c,p) == 1:
                    yield p,9,c,c


def check(executable, bound):
    result = run([str(Path(executable).resolve()),str(bound)], check=True,
                 text=True, capture_output=True)
    print(result.stderr.strip())
    rows = [tuple(map(int,line.split())) for line in result.stdout.splitlines()]
    assert len(rows) == len(set(rows))
    small = min(bound,5000)
    expected = {shape for shape in small_shapes(small) if accepted_inputs(*shape)}
    assert expected == {shape for shape in rows if shape[0] <= small}

    counts = Counter()
    for p,v,c,e in rows:
        inputs = accepted_inputs(p,v,c,e)
        assert inputs
        for q in inputs:
            for slopes in normalized_rows(p,v,c,e):
                x,y,z = (a*q+(i==2) for i,a in enumerate(slopes))
                assert (x | y | z) == x+y+z
            # The two full prime powers are distinct, and the fixed one can
            # be chosen as the smaller. Also enforce the exact assignment residues.
            residue, modulus = (587,1008) if v == 7 else (1267,1584)
            admissible = q >= p and (p*q-residue) % modulus == 0 and (c*q-1) % p == 0
            assert not admissible
            counts[v,q] += 1
    if bound == 500000:
        assert counts == {(7,1):40, (7,49):7, (9,1):1}
    print("Complete binary inputs by (multiplier, other factor):", sorted(counts.items()))
    print("All possible shapes directly rechecked through fixed power", small)


if __name__ == "__main__":
    cli = ArgumentParser(description=__doc__)
    cli.add_argument("executable")
    cli.add_argument("--bound", type=int, default=500000)
    args = cli.parse_args()
    check(args.executable,args.bound)
