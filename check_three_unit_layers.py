"""Verify three-unit predecessor layers and the finite table of possible grouped cuts."""

from collections import defaultdict
from itertools import combinations, combinations_with_replacement, permutations
from math import gcd, isqrt

from check_endpoint_pair_bound import prime_powers, rows
from check_seven_mixed import log_intervals, separating_denominator, valuation_cap


PRIMES = [p for p in range(2, 256) if all(p % k for k in range(2, isqrt(p) + 1))]
EXPONENT_BOUND = 100


def digit_units(n, base):
    result = []
    place = 1
    while n:
        n, digit = divmod(n, base)
        result.extend([place] * digit)
        place *= base
    return result


def free_part(value, primes):
    for p in primes:
        while value % p == 0:
            value //= p
    return value


def check_separation():
    logs = log_intervals(PRIMES)
    largest_denominator = largest_lifting_factor = 0
    for p, q in permutations(PRIMES, 2):
        denominator = separating_denominator(logs[p], logs[q])
        if q == 2:
            cap = (p + 1) & -(p + 1)
        else:
            cap = max(1, valuation_cap(p, q))
        largest_denominator = max(largest_denominator, denominator)
        largest_lifting_factor = max(largest_lifting_factor, cap)
    assert largest_lifting_factor == 28561
    assert 3 * largest_lifting_factor * EXPONENT_BOUND * 2**71 < 2**EXPONENT_BOUND
    return largest_denominator


def grouped_table():
    by_value = defaultdict(list)
    for p in PRIMES:
        powers = [p**i for i in range(EXPONENT_BOUND)]
        for a, b in combinations_with_replacement(powers, 2):
            w = a + b
            if (w + 1) % p and (a != b or p != 2):
                by_value[w].append(p)
    result = {(p, q, w) for w, bases in by_value.items() for p, q in combinations(bases, 2)
              if free_part(w, (p, q)) >= 3}

    # Independent membership test. Extract the two lowest-place valuations,
    # then test the remaining quotient for being one plus a prime power.
    def two_units(w, q):
        unit = free_part(w, (q,))
        if unit == 1:
            return False
        return free_part(unit - 1, (q,)) == 1 and (q != 2 or unit > 2)

    other = set()
    for p, q in combinations(PRIMES, 2):
        for i in range(EXPONENT_BOUND):
            for j in range(i if p > 2 else i + 1, EXPONENT_BOUND):
                w = p**i + p**j
                if gcd(w + 1, p*q) == 1 and two_units(w, q) and free_part(w, (p, q)) >= 3:
                    other.add((p, q, w))
    assert result == other
    assert len(result) == 25 and len({(p, q) for p, q, w in result}) == 23
    by_w = defaultdict(set)
    for p, q, w in result:
        by_w[w].add((p, q))
    triangles = {(p, q, r, w) for w, edges in by_w.items()
                 for p, q, r in combinations(sorted({v for pair in edges for v in pair}), 3)
                 if {(p, q), (p, r), (q, r)} <= edges}
    assert triangles == {(2, 11, 131, 132)}
    return sorted(result)



def orbit_exponents(table):
    """Upper bounds in thousandths for a block starting with any multiplication carry."""
    bounds = {}
    for smaller, larger, w in table:
        for p in (smaller, larger):
            if (p, w) in bounds:
                continue
            if w % p:
                assert (w + 1) % p and p > 2
                length, allowed = 1, (p + 1)//2
                for carry in range(w):
                    assert sum(d + (w*d+carry) % p < p for d in range(p)) == allowed
            else:
                transitions = [[(w*d+c)//p for d in range(p - c % p)] for c in range(w)]
                length = 64
                ways = [1]*w
                for level in range(1, length+1):
                    ways = [sum(ways[next_c] for next_c in choices) for choices in transitions]
                    if level == 2:
                        for c in (0, w//2, w-1):
                            # An independent whole-integer test at both place values.
                            direct = sum(all(n % m + (w*n+c) % m < m for m in (p, p*p))
                                         for n in range(p*p))
                            assert ways[c] == direct
                allowed = max(ways)
            low, high = 0, 1000
            target = allowed**1000
            while high - low > 1:
                mid = (low + high)//2
                if target <= p**(length*mid):
                    high = mid
                else:
                    low = mid
            bounds[p, w] = high
    pairs = [(bounds[p,w], bounds[q,w]) for p,q,w in table]
    result = max(max(pair) for pair in pairs), max(map(sum, pairs)), max(map(min, pairs))
    assert result == (876, 1619, 776)
    return result


def check_binomial_layers():
    # Reconstruct maximal layers directly from W. Use factorial valuations
    # to test every cut, independently of the affine comparison.
    for W, bases in [(12, (2, 3)), (105, (5, 7)), (4930, (2, 5, 17)),
                     (59052, (2, 3, 7))]:
        full = dict(prime_powers(W))
        layers = {}
        for p in bases:
            P = full[p]
            t = W//P-1
            assert t > 0 and t % p == 0
            layers[p] = P, free_part(t, (p,))
        edges = []
        for p, q in combinations(bases, 2):
            P, t = layers[p]
            Q, u = layers[q]
            if t*u*(P+Q+1) < W:
                edges.append((p, q))
        assert len(edges) >= len(bases)-1
        zeros = {}
        for p in bases:
            def factorial_valuation(x):
                total = 0
                while x:
                    x //= p
                    total += x
                return total
            vals = [factorial_valuation(x) for x in range(W+2)]
            zeros[p] = {a for a in range(2, W) if vals[a]+vals[W+1-a] == vals[W+1]}
        for p, q in edges:
            P, t = layers[p]
            Q, u = layers[q]
            difference = abs(P-Q)
            possible = set()
            if gcd(t, u) % difference == 0:
                a, remainder = divmod(W-min(P,Q), difference)
                assert remainder == 0
                possible = {a, W+1-a}
            assert zeros[p] & zeros[q] <= possible
        actual = set.intersection(*(zeros[p] for p in bases))
        assert actual == {12: {4, 9}, 105: {50, 56}}.get(W, set())


def exponent_class(g):
    t = 2**g
    d = 1 + t*t
    h = d + t
    k = 4*g + 1
    assert d*h <= 2**k and digit_units(h, 2) == [1, t, t*t]
    assert h > d and d <= min(5*h, (1+t)*(4*(1+t)+1))
    exponent = 0
    for bits in range(3, k):
        modulus = 2**(bits+1)
        if (d*pow(3, exponent, modulus) - h) % modulus:
            exponent += 2**(bits-2)
        assert (d*pow(3, exponent, modulus) - h) % modulus == 0
    period = 2**(k-2)
    assert pow(3, period, 2**k) == 1 and pow(3, period//2, 2**k) != 1
    assert exponent % 2 == 0 and 0 < exponent < period
    return exponent, period


def low_layer(value, p, d):
    # Once the first three units have appeared, the next nonzero digit can
    # only invalidate the unit budget. Allow trailing zeros above the number.
    modulus = 1
    while modulus < p * d * value:
        modulus *= p
        h = value % modulus
        units = digit_units(h, p)
        if len(units) > 3:
            break
        if d*h <= modulus and len(units) == 3:
            return h, units
    return None


def check_actual_coefficients(bound, possibilities):
    counts = [0, 0]
    weights = coefficient_count = 0
    for n in range(3, bound+1):
        powers = list(prime_powers(n-1))
        one = two = False
        for (p, P), (q, Q) in permutations(powers, 2):
            d = (n-1) // (P*Q)
            layer = low_layer(d*Q, p, d)
            if layer is None:
                continue
            h, (a, b, c) = layer
            ratio = min(b//a, c//b)
            one |= d % h != 0 and d <= min(5*h, (1+ratio)*(4*(1+ratio)+1))
            other = low_layer(d*P, q, d)
            if other:
                two |= h != other[0] or d % h != 0 or (min(p, q), max(p, q), h-1) not in possibilities
        counts[0] += one
        counts[1] += two
        if not (one or two):
            continue
        weights += 1
        # Every triple of actual coefficients has nontrivial gcd exactly when
        # every three inclusion-minimal prime supports meet.
        primes = [p for p in PRIMES if p <= n]
        minimal = set()
        for _, coefficient in rows(n):
            mask = sum(1 << i for i, p in enumerate(primes) if coefficient % p == 0)
            coefficient_count += 1
            if any(old & mask == old for old in minimal):
                continue
            minimal = {old for old in minimal if old & mask != mask}
            minimal.add(mask)
        assert all(a & b & c for a, b, c in combinations_with_replacement(minimal, 3))
    return counts, weights, coefficient_count


if __name__ == '__main__':
    print('Largest certified rational denominator:', check_separation())
    table = grouped_table()
    print('Necessary (p, q, s-1) triples:', table)
    print('Marginal, square, synchronized exponents in thousandths:', orbit_exponents(table))
    for g in (3, 4, 9, 27, 81):
        print('g, representative, period:', g, *exponent_class(g))
    assert [c for c in range(2048) if 65*pow(3, c, 8192) % 8192 == 73] == [386]
    check_binomial_layers()
    direct = check_actual_coefficients(250, set(table))
    assert direct == ([123, 23], 123, 217620)
    print('One-sided and paired weights, union, exact coefficients:', direct)
