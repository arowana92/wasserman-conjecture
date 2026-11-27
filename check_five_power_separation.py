"""Check exponent separation and finite ranges at cofactor five, k <= 35."""

from collections import defaultdict, deque
from functools import lru_cache
from itertools import combinations
from math import gcd, isqrt

from check_factor_five_binary_exponents import graph

# k, coefficient of Q in the middle part, required low k+1 bits of Q.
SHAPES = (
    (31, 4_192_014_075, 741_526_579),
    (32, 6_285_315, 5_727_322_795),
    (34, 34_312_377_675, 429_493_347),
)


def binary_graph(shape):
    k, c, q = shape
    starts, reverse, terminal = graph(k, c, q, 2, whole=True)
    assert len(starts) == 1
    # Index precisely the states that can reach a stopping state.
    queue = list(terminal)
    index = {s: i for i, s in enumerate(queue)}
    edges = [[] for s in queue]
    for t in queue:
        for s, digit in reverse.get(t, ()):
            if s not in index:
                index[s] = len(queue)
                queue.append(s)
                edges.append([])
            edges[index[s]].append((digit, index[t]))
    return index[starts[0]], edges, set(range(len(terminal)))


def prefixes(shape, data, bits):
    k, _, q = shape
    start, edges, _ = data
    words = [(start, q)]
    for j in range(k+1, bits):
        words = [(t, value+(digit << j)) for s, value in words for digit, t in edges[s]]
    return [value for _, value in words]


def log_three(r, bits):
    # The class is 3 mod 8. Adding 2^(j-2) to the exponent toggles bit j.
    modulus = 1 << bits
    power, exponent = 3, 1
    for j, step in enumerate(log_steps(bits), 3):
        if (r ^ power) & (1 << j):
            exponent += 1 << (j-2)
            power = power*step % modulus
    assert power == r % modulus and r % 8 == 3
    return exponent


@lru_cache(maxsize=None)
def log_steps(bits):
    return [pow(3, 1 << (j-2), 1 << bits) for j in range(3, bits)]


def shortest_vector(ratio, modulus):
    # Lattice of (a,b) with b = ratio*a mod modulus.
    a, b, c, d = 1, ratio, 0, modulus
    while True:
        u, v = a*a+b*b, c*c+d*d
        if v < u:
            a, b, c, d = c, d, a, b
            u = v
        dot = a*c+b*d
        if 2*abs(dot) <= u:
            return (a, b) if a > 0 else (-a, -b)
        quotient = (2*dot+u)//(2*u)
        c -= quotient*a
        d -= quotient*b


def separation(graphs):
    bound = 10**18
    coarse, fine = 124, 176
    n, m = 1 << (coarse-2), 1 << (fine-2)
    assert n > 2*bound**2
    by_prefix = defaultdict(set)
    counts = []
    for shape, data in zip(SHAPES, graphs):
        words = prefixes(shape, data, fine)
        logs = [log_three(q, fine) for q in words]
        assert len(set(logs)) == len(logs)
        residues = {x % n for x in logs}
        assert residues.isdisjoint(by_prefix)
        counts.append((len(residues), len(logs)))
        for x in logs:
            by_prefix[x % n].add(x)
    print("Logarithm classes at 124 and 176 binary digits:", counts)
    inverse = {x: pow(x, -1, n) for x in by_prefix}
    coarse_relations = 0
    for x, y in combinations(by_prefix, 2):
        a, b = shortest_vector(y*inverse[x] % n, n)
        if not (1 <= a <= bound and 1 <= b <= bound and a*b % 2 and gcd(a,b) == 1):
            continue
        assert (b*x-a*y) % n == 0
        coarse_relations += 1
        ratio = b*pow(a,-1,m) % m
        assert all(ratio*z % m not in by_prefix[y] for z in by_prefix[x])
    assert coarse_relations == 3947
    print("Short positive odd ratios at the two precisions:", coarse_relations, 0)


def path_bounds(graphs):
    counts = []
    for data, exponent in zip(graphs, (19,74,105)):
        _, edges, _ = data
        destinations = [[v for _,v in row] for row in edges]
        assert all(destinations)
        paths = [1]*len(edges)
        for _ in range(2048):
            paths = [sum(paths[v] for v in row) for row in destinations]
        maximum = max(paths)
        assert maximum**1000 < 2**(2048*exponent)
        counts.append(maximum)
    print("Maximum numbers of 2048-digit paths:", counts)


def primes_through(bound):
    sieve = bytearray([1])*(bound+1)
    sieve[:2] = b'\0\0'
    for p in range(2, isqrt(bound)+1):
        if sieve[p]:
            sieve[p*p::p] = bytes((bound-p*p)//p+1)
    return (n for n in range(2, bound+1) if sieve[n])


def avoids_layer_cut(total, a, b, p, digits):
    # The third part is total-a-b. Test successive factorial-floor layers.
    c = total-a-b
    power = 1
    for _ in range(digits):
        power *= p
        if total//power > a//power+b//power+c//power:
            return False
    return True


def prime_box(graphs):
    bits = 130
    modulus = 1 << (bits-2)
    prime_list = [r for r in primes_through(1_000_000) if r > 5 and r % 8 == 3]
    totals = []
    for shape, data in zip(SHAPES, graphs):
        k, c, _ = shape
        P = 1 << k
        logs = [log_three(x, bits) for x in prefixes(shape, data, bits)]
        bases = candidates = small = 0
        for r in prime_list:
            # At the base of Q=r^e, Q shifts all three slopes past the offset 1.
            if not avoids_layer_cut(5*P, 3*P, c, r, (5*P).bit_length()):
                continue
            bases += 1
            inv = pow(log_three(r, bits), -1, modulus)
            for x in logs:
                e = x*inv % modulus
                assert e % 2
                candidates += 1
                residue = pow(r, e % 12, 45)
                if P*residue % 15 != 7 or (5*P*residue+1) % 9 == 0:
                    continue
                small += 1
                for p, depth in ((3,12), (5,12), (3,32), (5,32)):
                    limit = p**depth
                    q = pow(r, e % ((p-1)*(limit//p)), limit)
                    if not avoids_layer_cut(1+5*P*q, 3*P*q, c*q, p, depth):
                        break
                else:
                    raise AssertionError((k,r,e))
        totals.append((bases,candidates,small))
    print("Prime bases, exponents and survivors modulo 45 before carry tests:", totals)


def finite_scalars(graphs):
    # Split each binary word in two. At fixed length high zero padding is unique.
    totals = []
    for shape, data, bits in zip(SHAPES, graphs, (1152,384,384)):
        k, c, q0 = shape
        P = 1 << k
        start, edges, terminal = data
        cut = (bits+k+1)//2
        low = [(start,q0)]
        for j in range(k+1,cut):
            low = [(t, a+(b << j)) for s,a in low for b,t in edges[s]]
        # Cache high suffixes only for states and lengths that the low halves use.
        @lru_cache(maxsize=None)
        def suffixes(s, length):
            if length == 0:
                return (0,) if s in terminal else ()
            return tuple(2*x+b for b,t in edges[s] for x in suffixes(t,length-1))

        depth = 7
        modulus = 3**depth*125
        ok = bytearray(modulus)
        for q in range(modulus):
            if P*q % 15 == 7 and (5*P*q+1) % 9:
                ok[q] = (avoids_layer_cut(1+5*P*q,3*P*q,c*q,3,depth)
                         and avoids_layer_cut(1+5*P*q,3*P*q,c*q,5,3))
        tested = passed = 0
        shift = pow(2, cut, modulus)
        @lru_cache(maxsize=None)
        def high(s):
            return [(v, v*shift % modulus) for v in suffixes(s,bits-cut)]

        for s, value in low:
            lower = value % modulus
            tails = high(s)
            tested += len(tails)
            for v, residue in tails:
                if not ok[(lower+residue) % modulus]:
                    continue
                passed += 1
                q = value+(v << cut)
                n, a, b = 1+5*P*q, 3*P*q, c*q
                assert any(not avoids_layer_cut(n,a,b,p,n.bit_length()) for p in (3,5))
        totals.append((bits,tested,passed))
        suffixes.cache_clear()
        high.cache_clear()
    print("Scalar bit limit, accepted binary words and low odd-base survivors:", totals)


if __name__ == "__main__":
    # The longest suffix in the first shape has 560 binary digits.
    import sys
    sys.setrecursionlimit(5000)
    graphs = [binary_graph(shape) for shape in SHAPES]
    separation(graphs)
    path_bounds(graphs)
    prime_box(graphs)
    finite_scalars(graphs)
