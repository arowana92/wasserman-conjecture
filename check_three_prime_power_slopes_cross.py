"""Independent residue search. All arithmetic, including the floor tests, is integral."""
from collections import defaultdict
from functools import cache
from itertools import combinations
from math import gcd


def powers(p, limit):
    value = p
    while value <= limit:
        yield value
        value *= p


def floor_layers(parts, base, stop=None):
    total = sum(parts)
    if stop is None:
        stop = total
    for m in powers(base, stop):
        if total // m != sum(part // m for part in parts):
            return False
    return True


def candidate_parts(own, other, prime, r):
    # Determine the multiplier modulo a power of r from the low digits of its part.
    total = own * other
    modulus, values = 1, [0]
    while modulus * r <= own:
        digit = total // modulus % r
        values = [v + t * modulus for t in range(digit + 1) for v in values]
        modulus *= r
    reciprocal = pow(other, -1, modulus)
    for low in values:
        x0 = low * reciprocal % modulus
        for x in range(x0 or modulus, own, modulus):
            if x % prime == 0:
                continue
            slope = x * other
            high = slope // modulus
            if floor_layers((high, total // modulus - high), r):
                yield x


@cache
def prime_factors(n):
    out = []
    t = 2
    while t * t <= n:
        if n % t == 0:
            out.append(t)
            while n % t == 0:
                n //= t
        t = 3 if t == 2 else t + 2
    if n > 1:
        out.append(n)
    return tuple(out)


@cache
def fixed_power(r, period, prime):
    h = r**period - 1
    modulus = 1
    while h % prime == 0:
        h //= prime
        modulus *= prime
    return modulus


def grouped_parts(own, other, prime, r, period):
    fixed = gcd(own, r**period - 1)
    groups = defaultdict(list)
    for x in candidate_parts(own, other, prime, r):
        groups[x * other % fixed].append(x)
    return groups, fixed


def surviving_cell(p, q, r, P, Q, period):
    # Index the two ranges by the residues prescribed by the third exponent.
    X, px = grouped_parts(P, Q, p, r, period)
    if not X:
        return []
    Y, qy = grouped_parts(Q, P, q, r, period)
    total = P * Q
    residues = defaultdict(list)
    for e in range(period):
        residues[pow(r, -e, px), pow(r, -e, qy)].append(e)
    hits = []
    for (u, v), exponents in residues.items():
        ys = Y.get(v, ())
        if not ys:
            continue
        for x in X.get(u, ()):
            rp = pow(x * Q, -1, P)
            for y in ys:
                g = gcd(x, y)
                while g % r == 0:
                    g //= r
                if g < 3 or x * Q + y * P >= total:
                    continue
                if not floor_layers((y * rp, (Q-y) * rp), p, P):
                    continue
                rq = pow(y * P, -1, Q)
                if not floor_layers((x * rq, (P-x) * rq), q, Q):
                    continue
                slopes = (x * Q, y * P, total - x * Q - y * P)
                if not floor_layers(slopes, r):
                    continue
                for ell in prime_factors(g):
                    required = 4 if ell == 2 else ell
                    if g % required:
                        continue
                    for e in exponents:
                        m = fixed_power(r, period, ell)
                        if (total * pow(r, e, m) + 1) % gcd(required, m):
                            continue
                        for base in (p, q, ell):
                            m = fixed_power(r, period, base)
                            z = pow(r, e, m)
                            row = (slopes[0]*z, slopes[1]*z, slopes[2]*z+1)
                            if not floor_layers(row, base, m):
                                break
                        else:
                            hits.append((P, Q, r, x, y, ell, e))
    return hits


def run():
    import sys
    sets = [(s, 10**7, 12) for s in combinations((2, 3, 5, 7), 3)]
    sets += [((2, 3, r), 10**6, 60) for r in (11, 13)]
    for triple, cap, period in sets:
        if len(sys.argv) > 1 and ','.join(map(str, triple)) not in sys.argv[1:]:
            continue
        cells = 0
        for p, q in combinations(triple, 2):
            r, = set(triple) - {p, q}
            for P in powers(p, cap):
                for Q in powers(q, cap):
                    hits = surviving_cell(p, q, r, P, Q, period)
                    assert not hits, hits
                    cells += 1
        print(triple, cells, 'families checked', flush=True)


if __name__ == '__main__':
    run()
