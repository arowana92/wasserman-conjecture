"""Check the local witnesses, exponent class, separation and interval on the 37-ray."""

from collections import deque
from math import gcd, lcm

from carry_blocks import CarryGraph
from check_higher_successors import LAST
from prime_certificates import primes_through

START = 2172940
PERIOD = 2718912
# Each integer satisfies the original shape congruence and has a carry-free row locally.
WITNESSES = {2: (443822116471505, 8), 7: (1676335771086440963, 3),
             17: (27225360927, 3), 67: (11694083, 1)}


def finishing_graph(r, slopes, offsets):
    """Keep precisely the edges that can be followed by a finite accepting suffix."""
    graph = CarryGraph(r, slopes, offsets)
    parents = [[] for _ in graph.states]
    for i, edges in enumerate(graph.edges):
        for _, j in edges:
            parents[j].append(i)
    live = {i for i in range(len(parents)) if graph.accept >> i & 1}
    todo = deque(live)
    while todo:
        for i in parents[todo.popleft()]:
            if i not in live:
                live.add(i)
                todo.append(i)
    return [{digit: j for digit, j in edges if j in live} for edges in graph.edges]


def read_prefix(edges, r, scalar, depth):
    state = 0
    for _ in range(depth):
        scalar, digit = divmod(scalar, r)
        if digit not in edges[state]:
            return None
        state = edges[state][digit]
    return state


def low_classes(r, depth):
    w = r**depth
    L = 48
    while pow(37, L, w) != 1:
        L *= r
    scalar = pow(37, 28, w)
    step = pow(37, 48, w)
    first = r != 2
    edges = finishing_graph(r, LAST.first if first else LAST.other,
                            (0, 0, 1) if first else (0, 0, 0))
    result = []
    for a in range(28, L, 48):
        v = scalar if first else (2278*scalar + 1)*pow(7, -1, w) % w
        if (r != 7 or (2278*scalar + 1) % 49) and read_prefix(edges, r, v, depth) is not None:
            result.append(a)
        scalar = scalar*step % w
    return L, result


def no_carry(parts, r):
    parts = tuple(parts)
    w = r
    while w <= sum(parts):
        if sum(parts)//w != sum(k//w for k in parts):
            return False
        w *= r
    return True


def local_witnesses(depth):
    residue, modulus = START, PERIOD
    for r, (scalar, k) in WITNESSES.items():
        w = r**k
        assert scalar % LAST.modulus == LAST.residue and (2278*scalar+1) % 49
        assert scalar % w == pow(37, START, w)
        assert pow(37, PERIOD, w) == 1 and pow(37, PERIOD, r*w) != 1
        if r == 2:
            V = (2278*scalar+1)//7
            assert no_carry([f*V for f in LAST.other], 2)
        else:
            x, y, z = LAST.first
            assert no_carry((x*scalar, y*scalar, z*scalar+1), r)
        # Extend the exponent, retaining all preceding primes' precision.
        while w < r**depth:
            w *= r
            digit = next(i for i in range(r) if pow(37, residue+i*modulus, w) == scalar % w)
            residue += digit*modulus
            modulus *= r
    for r, (scalar, _) in WITNESSES.items():
        w = r**depth
        P = pow(37, residue, w)
        if r == 2:
            V = (2278*P+1)*pow(7, -1, w) % w
            assert V == ((2278*scalar+1)//7) % w
            parts = tuple(f*V for f in LAST.other)
        else:
            x, y, z = LAST.first
            parts = x*P, y*P, z*P+1
        for _ in range(depth):
            assert sum(part % w for part in parts) < w
            w //= r
    print(f"Simultaneous witness through {depth} digits in all four bases; "
          f"the constructed exponent has {len(str(residue))} decimal digits.")

    return residue


def auxiliary_primes(bound, witness_exponent):
    # n in a = START + PERIOD*n. Use zero for its components outside the four bases.
    n = (witness_exponent-START)//PERIOD
    precision = 1
    for r, (_, h) in WITNESSES.items():
        power = 1
        while power < bound:
            power *= r
        assert power*r**h <= r**100
        precision *= power
    checked = 0
    for q in primes_through(bound):
        if q in WITNESSES or q == 37:
            continue
        m = q-1
        inside = gcd(m, precision)
        outside = m//inside
        t = n % inside
        t += inside * (-t*pow(inside, -1, outside) % outside)
        if (2278*pow(37, (START+PERIOD*t) % m, q) + 1) % q == 0:
            raise ArithmeticError(("additional weight prime", q))
        checked += 1
    print(f"All {checked} auxiliary primes through {bound} are excluded from the weight.")


def seven_classes(depth, edges):
    """Lift via the next admissible input digit, using the inverse of its digit increment."""
    if depth < 3:
        raise ValueError("the seed class already specifies three digits")
    final = 7**depth
    levels = []
    w, period = 343, PERIOD
    while w < final:
        step = pow(37, period, final)
        delta = pow(37, START, 7) * ((step-1)//w) % 7
        levels.append((w, period, pow(delta, -1, 7), [pow(step, i, final) for i in range(7)]))
        w *= 7
        period *= 7
    initial = pow(37, START, final)
    start_state = read_prefix(edges, 7, initial, 3)
    assert start_state is not None

    def lift(k, a, v, state):
        if k == len(levels):
            yield a
            return
        w, stride, inverse, steps = levels[k]
        current = v//w % 7
        for digit, child in edges[state].items():
            i = (digit-current)*inverse % 7
            yield from lift(k+1, a+i*stride, v*steps[i] % final, child)

    return period, lift(0, START, initial, start_state)


def check_separation(depth=30):
    edges = finishing_graph(7, LAST.first, (0, 0, 1))
    period, result = seven_classes(depth, edges)
    classes = sorted(result)
    gap = min(b-a for a, b in zip(classes, classes[1:] + [classes[0]+period]))
    assert gap > 0 and classes[-1] < period
    # Verify the direction of the logarithmic lifts independently at five new digits.
    small_period, small_classes = seven_classes(8, edges)
    direct = [a for a in range(START, small_period, PERIOD)
              if read_prefix(edges, 7, pow(37, a, 7**8), 8) is not None]
    assert direct == sorted(small_classes) and direct
    # Test actual representatives below one full period at the required bases.
    # First tabulate 18 low binary digits. Raising 37 to 2**16 fixes them.
    w = 2**18
    small = 2**16
    assert pow(37, small, w) == 1
    valid = set()
    for a in range(START % gcd(PERIOD, small), small, gcd(PERIOD, small)):
        V = (2278*pow(37, a, w)+1)*pow(7, -1, w) % w
        if (V & (2*V) | V & (4*V) | (2*V) & (4*V)) % w == 0:
            valid.add(a)
    count = 0
    for a in classes:
        if a % small not in valid:
            continue
        count += 1
        for r in (2, 17, 7, 67):
            w = r**100
            P = pow(37, a, w)
            if r == 2:
                V = (2278*P+1)*pow(7, -1, w) % w
                parts = tuple(f*V for f in LAST.other)
            else:
                x, y, z = LAST.first
                parts = x*P, y*P, z*P+1
            while w > 1 and sum(f % w for f in parts) < w:
                w //= r
            if w > 1:
                break
        else:
            raise ArithmeticError(("unresolved actual exponent", a))
    if depth == 30:
        assert (period, len(classes), gap, count) == (
            178666130578561716682765137216, 878080, 94438145595691968, 70603)
    print(f"Depth {depth}: {len(classes)} classes modulo {period}, minimum cyclic gap {gap}.")
    print(f"{count} representatives passed 18 binary digits; every exponent through {period} is excluded.")


if __name__ == "__main__":
    classes = [low_classes(r, k) for r, k in ((2, 9), (7, 3), (17, 3))]
    assert classes == [(384, [76, 268]), (2352, [2044]), (13872, [8908])]
    values = [(76, 192), (2044, 2352), (8908, 13872)]
    assert PERIOD == lcm(*(m for a, m in values)) and 0 <= START < PERIOD
    assert all(START % m == a for a, m in values)
    print(f"Every exceptional exponent equals {START} modulo {PERIOD}.")
    witness_exponent = local_witnesses(100)
    auxiliary_primes(10**7, witness_exponent)
    check_separation()
