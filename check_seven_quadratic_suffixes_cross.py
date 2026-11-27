"""Replay the residual input search on the CRT progression using squared endpoints."""

from argparse import ArgumentParser
from math import isqrt, lcm


def independent_graph(D, z, w):
    A = (2*z, 2*w, D-2*z-2*w)
    M = lcm(*A)
    cuts = {k for a in A for k in range(0, M, M//a)}
    possible = {(A[0]*k//M, A[1]*k//M, A[2]*k//M-c)
                for k in cuts for c in range(2)}

    def binary_ok(row):
        return min(row) >= 0 and sum(bin(v).count('1') for v in row) == bin(sum(row)).count('1')

    start = frozenset(s for s in possible if binary_ok(tuple(a+x for a,x in zip(A,s))))
    nodes, ids, transitions = [start], {start:0}, []
    columns = ((0,0,0),(1,0,0),(0,1,0),(0,0,1))
    for group in nodes:
        edges = []
        for digit in (0,1):
            # Invert the division of each part by two for each possible output digit.
            old = frozenset(p for s in group for c in columns
                            if (p := tuple(2*x+y-digit*a for x,y,a in zip(s,c,A))) in possible)
            if old:
                if old not in ids:
                    ids[old] = len(nodes)
                    nodes.append(old)
                edges.append((digit,ids[old]))
        transitions.append(edges)
    count = {0:1}
    seen = {}
    while (key := tuple(sorted(count.items()))) not in seen:
        seen[key] = sum(count.values())
        next_count = {}
        for node,mult in count.items():
            for _,dst in transitions[node]:
                next_count[dst] = next_count.get(dst,0)+mult
        count = next_count
    max_count = max(seen.values())
    assert max_count == (5 if D == 7383 else 1)
    return transitions


def input_grid(D, a, b, z, w):
    # First take P in its class modulo D. Intersect the two exact residue-one
    # constraints on R by searching only a steps and then b steps.
    P = next(p for p in range(D) if (10*p+7)%D == 0)
    step = D
    for group,v in ((a,z),(b,w)):
        for _ in range(group):
            R = (10*P+7)//D
            if (2*P+1)*R*v % group == 1:
                break
            P += step
        else:
            raise AssertionError('inconsistent bin')
        step *= group
    return 10*P+6, 10*step


def run(shape, bits):
    D,a,b,z,w = shape
    graph = independent_graph(D,z,w)
    base, stride = input_grid(*shape)
    start = 26
    low = 1 << start
    K = 160*D*low
    d = D.bit_length()-1
    words = [(0,1)]
    grid = []
    for length in range(1, start-d+3):
        if length >= start-d-1:
            for node,n in words:
                val = base + stride*((isqrt(K*n)-base)//stride)
                grid.append((node,val,K*n-val*val))
        words = [(target,2*n+digit) for node,n in words for digit,target in graph[node]]
    covered = odd = 0
    for j in range(start,bits):
        next_grid = []
        for node,val,delta in grid:
            step_square = stride*(2*val+stride)
            assert (val < 0 or delta >= 0) and delta < step_square and val+stride > 0
            if step_square <= delta+K:
                P = (val+stride-6)//10
                if max(10**8+1,low) <= P < 2*low:
                    covered += 1
                    F = (2*P+1)*((10*P+7)//D)
                    N = (2*P+1)*(5*P+1)
                    parts = (F*z, F*w, N-F*(z+w))
                    assert min(parts)>0
                    odd += sum(bin(i).count('1') for i in parts) == bin(N).count('1')
            # Squared distance from the next threshold to the point whose index
            # in the progression doubles. Advance it to the last point below.
            y = 2*val-base
            initial = 4*delta+base*(4*val-base)
            for digit,target in graph[node]:
                yy, ee = y, initial + 2*K*digit
                gap = stride*(2*yy+stride)
                while yy+stride <= 0 or ee >= gap:
                    ee -= gap
                    yy += stride
                    gap += 2*stride**2
                next_grid.append((target,yy,ee))
        grid = next_grid
        K *= 2
        low *= 2
    print(f'D={D}: independent residue progression found {covered} eligible inputs below 2^{bits}, '
          f'{odd} with a carry-free row.', flush=True)
    assert odd == 0
    return covered


if __name__ == '__main__':
    cli = ArgumentParser(description=__doc__)
    cli.add_argument('bits', type=int, default=65536, nargs='?')
    nbits = cli.parse_args().bits
    assert nbits >= 27
    for row in [(4587,11,417,1668,154),(7383,23,107,642,276),(9453,23,411,822,92)]:
        run(row, nbits)
