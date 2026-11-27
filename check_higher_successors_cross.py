"""Independently recover the first-row scale residues at all 23 final assignments."""

from collections import Counter
from functools import lru_cache
from math import gcd, lcm, prod

from flint import fmpz

from check_higher_successor_ray import PERIOD, START, WITNESSES, finishing_graph, no_carry, read_prefix
from check_higher_successors import LAST, assignments, extendable, local_cases
from check_successor_exponents import required_bases
from search_endpoint import atoms, smallest_factors


def scale_orbit(p, m):
    orbit = [1]
    v = p % m
    while v != 1:
        orbit.append(v)
        v = v*p % m
    return orbit


def compare_assignment(Q, p, height, u, T):
    d = Q-height
    factors = atoms(Q, smallest_factors(Q))
    B = Q//T
    while B % p == 0:
        B //= p
    m = B*u
    orbit = [v for v in scale_orbit(p, m) if (Q*v+1) % u == 0]
    all_scales = (1 << len(orbit)) - 1
    masks = {}
    for r, w in factors:
        if B % r:
            continue
        by_value = {}
        for i, v in enumerate(orbit):
            t = v % w
            by_value[t] = by_value.get(t, 0) | (1 << i)
        masks[w] = by_value

    @lru_cache(None)
    def permitted(w, x, y):
        z = (Q-x-y) % w
        bits = 0
        for t, block in masks[w].items():
            if sorted((t*x % w, t*y % w, (t*z+1) % w)) == [0, 0, 1]:
                bits |= block
        return bits

    direct = {}
    for y in range(u, d, u):
        for x in range(1, d-y+1):
            large = Q-x-y
            choices = [(y, large, x)]
            if y <= x:
                choices.append((y, x, large))
            for first in choices:
                a, b, c = first
                if b % u or not no_carry(first, p):
                    continue
                bits = all_scales
                for w in masks:
                    bits &= permitted(w, a % w, b % w)
                    if not bits:
                        break
                if bits:
                    direct[first] = bits
    by_factors = {}
    for v, S, _, first, R, M, _, _ in assignments(Q, p, height, factors):
        if (v, S) == (u, T):
            bits = sum(1 << i for i, value in enumerate(orbit) if value % M == R)
            by_factors[first] = by_factors.get(first, 0) | bits
    if direct != by_factors or not direct:
        raise ArithmeticError(("exhaustive scalar residues", Q, p, u, T))
    return len(direct)


@lru_cache(maxsize=4)
def graph_for(r, slopes, first):
    return finishing_graph(r, slopes, (0, 0, 1) if first else (0, 0, 0))


def independent_extension(c, r):
    first = r in required_bases(c)[0]
    graph = graph_for(r, c.first if first else c.other, first)
    w, depth = r, 1
    while w*r <= 4096:
        w *= r
        depth += 1
    period = scale_orbit(c.p, c.modulus)
    v = pow(c.p, period.index(c.residue), w)
    initial = v
    step = pow(c.p, len(period), w)
    while True:
        power = gcd(c.u, w)
        unitary = c.u % r or power == w or (c.Q*v+1) % (power*r)
        scalar = v if first else (c.Q*v+1)*pow(c.u, -1, w) % w
        if unitary and read_prefix(graph, r, scalar, depth) is not None:
            return True
        v = v*step % w
        if v == initial:
            return False


def check_auxiliary_primes():
    # Explicit discrete logarithms for n in a = START + PERIOD*n.
    conditions = {2: (16777216, 16730577), 7: (40353607, 2128588),
                  17: (24137569, 12821839), 67: (20151121, 297731)}
    for r, (m, n) in conditions.items():
        scalar, h = WITNESSES[r]
        assert m//r < 10**7 <= m and m % r == 0
        assert pow(37, START+PERIOD*n, m*r**h) == scalar % (m*r**h)
    small_primes = [2]
    outside_period = 1
    total = 0
    for q in range(3, 10**7+1, 2):
        if not fmpz(q).is_prime():
            continue
        if q <= 5000:
            small_primes.append(q)
        if q in conditions or q == 37:
            continue
        outside = q-1
        t = 0
        # Add CRT components separately, each divisible by all the other components.
        for r, (m, n) in conditions.items():
            power = gcd(q-1, m)
            outside //= power
            if power > 1:
                other = (q-1)//power
                t += other * (n*pow(other, -1, power) % power)
        if pow(37, START+PERIOD*(t % (q-1)), q)*2278 % q == q-1:
            raise ArithmeticError(("auxiliary-prime certificate", q))
        if q <= 5000:
            outside_period = lcm(outside_period, outside)
        total += 1
    assert total == 664574
    # At 5,000 also build a single exponent and test all primes by one gcd.
    n = 0
    stride = outside_period
    for m, value in conditions.values():
        n += stride*((value-n)*pow(stride, -1, m) % m)
        stride *= m
    a = START+PERIOD*n
    primorial = prod(small_primes)
    assert gcd(primorial, 2278*pow(37, a, primorial)+1) == 7
    print("Independent prime and CRT checks agree through ten million.")
    print("A single primorial gcd confirms all primes through 5,000.")


def solve():
    final = list(local_cases(1000, Counter()))
    keys = {(c.Q, c.p, c.height, c.u, c.T) for c in final}
    rows = sum(compare_assignment(*key) for key in sorted(keys))
    if (len(final), len(keys)) != (32, 23):
        raise ArithmeticError("unexpected classification")
    print(f"All scalar-residue sets agree for {rows} first rows at {len(keys)} assignments.", flush=True)
    by_base = Counter()
    others = []
    for c in final:
        first, second = required_bases(c)
        r = next((r for r in sorted(set(first+second)) if not extendable(c, r)), None)
        if r is None:
            others.append(c)
        elif independent_extension(c, r):
            raise ArithmeticError(("reverse graph admits a rejected prefix", c, r))
        else:
            by_base[r] += 1
    assert by_base == {2: 18, 3: 1, 7: 1, 41: 2, 83: 6, 101: 1}
    assert {c.Q for c in others} == {2278, 5447, 38141} and LAST in others
    for c in others:
        if c == LAST:
            assert all(independent_extension(c, r) for r in (2, 7, 17, 67))
            continue
        r = 2 if c.Q == 5447 else 3
        assert c.p % r == 1 and (c.Q+1) % r == 0 and c.u % r != 0
        assert all(sum(k*i % r for k in c.other) >= r for i in range(1, r))
    print("Reverse reachability confirms all 29 exclusions; fixed primes confirm two more.")


if __name__ == "__main__":
    solve()
    check_auxiliary_primes()
