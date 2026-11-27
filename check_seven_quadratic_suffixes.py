"""Count all binary suffixes of the three residual shapes and search their inputs."""

from argparse import ArgumentParser
from collections import defaultdict
from fractions import Fraction
from math import isqrt


RESIDUAL = ((4587,11,417,1668,154), (7383,23,107,642,276), (9453,23,411,822,92))


def upper_graph(D, z, w):
    slopes = (2*z, 2*w, D-2*z-2*w)
    breaks = {Fraction(i,a) for a in slopes for i in range(a)}
    offsets = {tuple(int(a*f) for a in slopes) for f in breaks}
    states = sorted({(x,y,s-e) for x,y,s in offsets for e in (0,1)})
    index = {s:i for i,s in enumerate(states)}
    back = [defaultdict(list), defaultdict(list)]
    terminal = set()
    for i, triple in enumerate(states):
        x,y,s = triple
        if s >= 0 and not ((x&y) | (x&s) | (y&s)):
            terminal.add(i)
        for bit in (0,1):
            values = [u + bit*a for u,a in zip(triple, slopes)]
            if sum(v%2 for v in values) <= 1:
                target = tuple(v//2 for v in values)
                back[bit][index[target]].append(i)

    def predecessor(group, bit):
        return frozenset(i for target in group for i in back[bit][target])

    first = predecessor(terminal, 1)
    assert first
    groups, known, graph = [first], {first:0}, []
    for group in groups:
        edges = []
        for bit in (0,1):
            source = predecessor(group, bit)
            if not source:
                continue
            if source not in known:
                known[source] = len(groups)
                groups.append(source)
            edges.append((bit, known[source]))
        graph.append(edges)
    return len(states), graph


def profiles(graph):
    profile = (1,) + (0,)*(len(graph)-1)
    seen, counts = {}, []
    while profile not in seen:
        seen[profile] = len(counts)
        counts.append(sum(profile))
        nxt = [0]*len(graph)
        for source, multiplicity in enumerate(profile):
            for _, target in graph[source]:
                nxt[target] += multiplicity
        profile = tuple(nxt)
    start = seen[profile]
    return start, len(seen)-start, max(counts), set(counts[start:]), counts


def input_class(D, a, b, z, w):
    R, step = 7*pow(D,-1,10)%10, 10
    for group, part in ((a,z),(b,w)):
        residue = -5*pow(2*part,-1,group)%group
        R += step*((residue-R)*pow(step,-1,group)%group)
        step *= group
    return (D*R-7)//10, D*(step//10)


def start_words(graph, lengths, K):
    words, result = [(1,0)], []
    for length in range(1, max(lengths)+1):
        if length in lengths:
            for n, node in words:
                s = isqrt(K*n)
                result.append((node,n,s,K*n-s*s))
        words = [(2*n+bit,target) for n,node in words for bit,target in graph[node]]
    return result


def advance(graph, words, K):
    result = []
    for node, n, s, error in words:
        for bit, target in graph[node]:
            # 2K(2n+bit) = (2s)^2 + 4*error + 2K*bit.
            extra = 4*error + 2*K*bit
            correction = extra//(4*s)
            remaining = extra - correction*(4*s+correction)
            while remaining < 0:
                remaining += 4*s + 2*correction - 1
                correction -= 1
            root = 2*s + correction
            assert 0 <= remaining <= 2*root
            result.append((target,2*n+bit,root,remaining))
    return result


def high_range(shape, graph, bits):
    D, a, b, z, w = shape
    residue, step = input_class(*shape)
    start = (10**8).bit_length()-1
    offset = D.bit_length()-1
    floor = 1 << start
    K = 160*D*floor
    words = start_words(graph, range(start-offset-1, start-offset+3), K)
    tested = in_class = odd = 0
    for j in range(start, bits):
        low, high = max(10**8+1, floor), 2*floor-1
        for _, n, s, error in words:
            tested += 1
            # (10P+6)^2 lies in (Kn,K(n+1)].
            first = max(low, (s-6)//10+1)
            P = first + (residue-first)%step
            if P > high:
                continue
            gap = 10*P + 6 - s
            if gap*(2*s+gap) > K+error:
                continue
            in_class += 1
            T, R = 2*P+1, (10*P+7)//D
            x, y = T*R*z, T*R*w
            c = T*(5*P+1) - x-y
            assert c > 0 and (10*P+6)**2-1 >= K*n
            if not ((x&y) | (x&c) | (y&c)):
                odd += 1
        # Recompute the square-root invariant directly at spaced checkpoints.
        if j == start or j+1 == bits or j%8192 == 0:
            for _, n, s, error in words:
                assert K*n == s*s+error and s == isqrt(K*n)
        words = advance(graph, words, K)
        floor *= 2
        K *= 2
    print(f'D={D}: {tested} upper words, {in_class} integer inputs in the bin congruences, '
          f'{odd} odd companions below 2^{bits}.', flush=True)
    assert odd == 0


def check_other_powers():
    # At the first bin each residual congruence makes R a nonsquare.
    for D,a,b,z,w in RESIDUAL:
        assert not any((2*z*r*r+5)%a == 0 for r in range(a))
        if D%10 == 7:
            continue
        count = above = 0
        for q in range(9,isqrt(7*D)+1,10):
            if any(q%p == 0 for p in range(2,isqrt(q)+1)):
                continue
            e = 3
            while q**(e-1) <= 7*D:
                count += 1
                R = q**e
                if (D*R-7)//10 > 10**8:
                    above += 1
                    assert (2*R*z+5)%a or (2*R*w+5)%b
                e += 2
        print(f'D={D}: {count} bounded odd powers, {above} above the finite input theorem.')


if __name__ == '__main__':
    cli = ArgumentParser(description=__doc__)
    cli.add_argument('bits', type=int, default=65536, nargs='?')
    args = cli.parse_args()
    assert args.bits >= 27
    check_other_powers()
    expected = {4587:(8222,138,0,138,1,1), 7383:(9192,283,28,1166,5,4),
                9453:(11272,68,0,68,1,1)}
    for shape in RESIDUAL:
        D,a,b,z,w = shape
        size, graph = upper_graph(D,z,w)
        pre, period, top, repeating, counts = profiles(graph)
        assert (size,len(graph),pre,period,top,next(iter(repeating))) == expected[D]
        assert len(repeating) == 1
        if D == 7383:
            assert [j+1 for j,v in enumerate(counts) if v>4] == [28]
        print(f'D={D}: {size} offsets, {len(graph)} subsets, {pre} initial profiles, '
              f'period {period}, largest count {top}.', flush=True)
        high_range(shape, graph, args.bits)
