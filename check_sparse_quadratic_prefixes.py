"""Restrict the leading digits and power exponents of the twelve quadratic rows."""
from fractions import Fraction
from math import gcd, isqrt

from check_sparse_quadratic_powers import KNOWN_SCALE, RESIDUAL_ROWS
from prime_certificates import primes_through


def carry_graph(d, i, j):
    slopes = (3*i, 3*j, d-3*(i+j))
    states = [None]  # The root must read a nonzero digit of the total.
    index = {}
    edges = []
    for state in states:
        z = state or (0, 0, 0)
        following = []
        for position in range(1 if state is None else 0, 4):
            digit = int(position != 0)
            target = tuple(2*z[k]+d*(position == k+1)-s*digit
                           for k, s in enumerate(slopes))
            if all(s-d <= v <= s for s, v in zip(slopes, target)):
                if target not in index:
                    index[target] = len(states)
                    states.append(target)
                following.append((index[target], digit))
        edges.append(following)
    return edges


def limiting_totals(edges):
    # Repeatedly delete dead ends. Every walk that leaves the remaining set
    # terminates within the number of vertices.
    live = set(range(len(edges)))
    while True:
        new = {a for a in live if any(b in live for b, _ in edges[a])}
        if new == live:
            break
        live = new
    graph = {a: [(b, bit) for b, bit in edges[a] if b in live] for a in live}
    limits = set()
    path = []
    levels = {}

    def visit(a, value):
        if a in levels:
            start = levels[a]
            # A cycle has a unique live continuation, even at later visits.
            assert all(len(graph[v]) == 1 for v in path[start:])
            width = len(path)-start
            prefix, block = divmod(value, 1 << width)
            limits.add((prefix+Fraction(block, (1 << width)-1))/(1 << start))
            return
        levels[a] = len(path)
        path.append(a)
        for b, bit in graph[a]:
            visit(b, 2*value+bit)
        path.pop()
        del levels[a]

    visit(0, 0)
    assert limits and all(Fraction(1, 2) <= f <= 1 for f in limits)
    return limits


def interval_check(row, lower_p, upper_p):
    d, i, j = row
    edges = carry_graph(*row)
    frontier = {(0, 0)}
    qmin = max(2, (3*lower_p-1+d-1)//d)
    qmax = (3*upper_p-1)//d
    scales = integral = 0
    odd = set()
    for h in range(1, qmax.bit_length()):
        frontier = {(b, 2*value+bit) for a, value in frontier for b, bit in edges[a]}
        lo = max(1 << h, qmin)
        hi = min((2 << h)-1, qmax)
        if lo > hi:
            continue
        left = ((d*lo+2)**2-1)//3
        right = ((d*hi+2)**2-1)//3
        for n in range(left.bit_length(), right.bit_length()+1):
            for prefix in {value for _, value in frontier}:
                scale = 1 << (n-h)
                # (d*q+2)^2 = 3*N+1 lies strictly between the lower
                # endpoint 3*prefix*scale and the upper endpoint minus one.
                first = (isqrt(3*prefix*scale)+d-2)//d
                last = (isqrt(3*(prefix+1)*scale-2)-2)//d
                for q in range(max(lo, first), min(hi, last)+1):
                    scales += 1
                    p, rem = divmod(d*q+1, 3)
                    if rem or q % 2 == 0:
                        continue
                    integral += 1
                    v = d*q+3
                    a, b = v*q*i, v*q*j
                    if a & b == 0 and (a | b) & ~(p*v) == 0:
                        odd.add(q)
    return scales, integral, odd


def controls():
    count = 0
    for row in ((2044,128,512), (2626,1,4), (2984,64,257)):
        d, i, j = row
        hi = (d*1000+1)//3
        _, _, found = interval_check(row, 1, hi)
        direct = set()
        for q in range(3, 1001, 2):
            p, rem = divmod(d*q+1, 3)
            if rem or p > hi:
                continue
            v = d*q+3
            parts = (q*i*v, q*j*v, (p-q*(i+j))*v)
            if sum(bin(x).count("1") for x in parts) == bin(p*v).count("1"):
                direct.add(q)
        assert found == direct and found
        count += len(found)
    print("Direct binary controls recovered:", count)


def homogeneous(c, base):
    for value in (c.numerator, c.denominator):
        for r in (2, base):
            while value % r == 0:
                value //= r
        if value != 1:
            return False
    return True


def exponent_check(constants):
    from flint import arb, ctx
    ctx.prec = 600
    log2 = arb(2).log()
    shifts = {c: (arb(c.numerator)/c.denominator).log()/log2 for c in constants}

    def nearest(x):
        return int((x+arb(1)/2).floor().unique_fmpz())

    base_count = resonances = max_bits = 0
    for base in primes_through(1_000_000):
        if base == 2:
            continue
        base_count += 1
        alpha = arb(base).log()/log2
        den = 1 << 490
        num = nearest(den*alpha)
        # Convergents propose denominators. Ball arithmetic checks every
        # error bound independently of the choice of rational approximation.
        up, u, tp, t = 0, 1, 1, 0
        pending = list(constants)
        while den and pending:
            a, rem = divmod(num, den)
            num, den = den, rem
            up, u = u, a*u+up
            tp, t = t, a*t+tp
            if t < 1 << 83:
                continue
            assert t < 1 << 360 and gcd(t, u) == 1
            if not (abs(t*alpha-u)*(1 << 78) < 1):
                continue
            max_bits = max(max_bits, t.bit_length())
            remaining = []
            for c in pending:
                if homogeneous(c, base):
                    resonances += 1
                    continue
                value = t*shifts[c]
                if not (abs(value-nearest(value)) > arb(1)/16):
                    remaining.append(c)
            pending = remaining
        assert not pending, (base, pending)
    assert (base_count, resonances, max_bits) == (78497, 4, 117)
    print("Odd prime bases:", base_count, "homogeneous cases:", resonances,
          "largest denominator bits:", max_bits)


if __name__ == "__main__":
    controls()
    constants = set()
    scales = integral = max_states = 0
    for row in RESIDUAL_ROWS:
        d, i, j = row
        edges = carry_graph(*row)
        max_states = max(max_states, len(edges))
        for f in limiting_totals(edges):
            constants.update((3/f, d*d/(3*f)))
        n, k, odd = interval_check(row, KNOWN_SCALE+1, 10**1000)
        assert not odd
        scales += n
        integral += k
    assert max_states <= 305 and len(constants) == 43
    assert max(max(c.numerator, c.denominator).bit_length() for c in constants) == 42
    print("Maximum graph size:", max_states, "rational constants:", len(constants))
    print("Prefix candidates:", scales, "integral odd scales:", integral, "odd rows: 0")
    assert (scales, integral) == (156259, 25995)
    exponent_check(constants)
