"""Classify the phases of four prime radicals and check exceptional gcd divisors."""
from itertools import permutations
from math import gcd, lcm
from random import Random

from flint import fmpz_poly
from sympy import isprime

from check_binomial_sparse_gcd import X, primitive


# Heights are doubled, sorted and translated to start at zero. For each
# height pattern the last number counts ordered minimal relations with
# the phase of the first term fixed at one.
CASES = {
    2: ((840, 9999999241, 140347707, 24),
        {(0, 0, 0, 1): 12, (0, 0, 0, 2): 6, (0, 0, 0, 3): 0,
         (0, 0, 1, 1): 6, (0, 0, 2, 2): 2, (0, 0, 1, 2): 4,
         (0, 0, 1, 3): 2, (0, 0, 2, 3): 2, (0, 0, 2, 4): 1}),
    3: ((420, 9999999661, 407226638, 12),
        {(0, 0, 0, 1): 6, (0, 0, 0, 2): 1, (0, 0, 1, 1): 4, (0, 0, 1, 2): 4}),
}


def phase_controls(p, parameters, cases):
    order, prime, root, final_order = parameters
    assert isprime(prime) and pow(root, order, prime) == 1
    assert all(pow(root, order//q, prime) != 1 for q in (2, 3, 5, 7))
    powers = [pow(root, j, prime) for j in range(order)]
    radical_period = 8 if p == 2 else 12
    shift = order//radical_period
    radical = (powers[shift]+powers[-shift]) % prime
    assert radical**2 % prime == p
    # The lowest two terms are 1 and a phase. Their nonzero sum
    # determines the only possible second phase for any upper pair.
    lower = {(-1-x) % prime: j for j, x in enumerate(powers) if j != order//2}
    phi = fmpz_poly.cyclotomic(order)
    phase = [(X**i) % phi for i in range(order)]
    rho = phase[shift]+phase[-shift]
    assert rho**2 % phi == p
    verified = 0
    exact_orders = set()
    for heights, expected in cases.items():
        h, k = heights[2:]
        upper = [pow(radical, k, prime)*x % prime for x in powers]
        first = [pow(radical, h, prime)*x % prime for x in powers]
        count = 0
        for a, x in enumerate(first):
            for b, y in enumerate(upper):
                j = lower.get((x+y) % prime)
                if j is None:
                    continue
                indices = (0, j, a, b)
                if any(heights[i] == heights[t] and (indices[i]-indices[t]) % order == order//2
                       for i in range(4) for t in range(i)):
                    continue
                if (1+phase[j]+rho**h*phase[a]+rho**k*phase[b]) % phi:
                    continue
                actual_order = order//gcd(order, j, a, b)
                assert final_order % actual_order == 0
                exact_orders.add(actual_order)
                count += 1
        assert count == expected, (p, heights, count, expected)
        verified += count
    assert lcm(*exact_orders) == final_order
    print(f"Base {p}: {verified} ordered relations, phase orders {sorted(exact_orders)}")



def direct_phase_counts(p, period, cases):
    # Enumerate directly in characteristic zero at the proved final order.
    # The modulus sieve above is still needed to exclude larger orders.
    phi = fmpz_poly.cyclotomic(period)
    phase = [(X**j) % phi for j in range(period)]
    shift = period//(8 if p == 2 else 12)
    radical = phase[shift]+phase[-shift]
    for heights, expected in cases.items():
        admissible = 0
        seen_orders = set()
        for j, second in enumerate(phase):
            for a, third in enumerate(phase):
                target = -(1+second+radical**heights[2]*third) % phi
                for b, fourth in enumerate(phase):
                    if target != (radical**heights[3]*fourth) % phi:
                        continue
                    indices = (0, j, a, b)
                    if any(heights[i] == heights[t] and (indices[i]-indices[t]) % period == period//2
                           for i in range(4) for t in range(i)):
                        continue
                    admissible += 1
                    seen_orders.add(period//gcd(period, j, a, b))
        assert admissible == expected
        yield heights, sorted(seen_orders)

def offsets(cases):
    return {tuple(v-a[0] for v in a[1:])
            for shape, count in cases.items() if count
            for a in permutations(shape)}


def integer_controls(p, period, candidates):
    rng = Random(785+p)
    count = strict = 0
    for _ in range(1300):
        d = rng.randrange(1, 40)
        positions = (0, *sorted(rng.randrange(50) for _ in range(3)))
        sign = rng.choice((-1, 1))
        original = (X**d+sign).gcd(sum(X**j for j in positions))
        base = abs(int(original(p)))
        actual = gcd(p**d+sign, sum(p**j for j in positions))
        outside = actual//base
        assert actual == base*outside
        allowed = base
        for delta in candidates:
            scale = gcd(2*d, *(2*v-t for v, t in zip(positions[1:], delta)))
            allowed = lcm(allowed, base*gcd(outside, p**(period*scale//2)-1))
        for slopes in ((d, *positions[1:]), tuple(rng.randrange(12) for _ in range(4))):
            f = p**d*X**slopes[0]+sign
            g = 1+sum(p**j*X**v for j, v in zip(positions[1:], slopes[1:]))
            line = abs(int(primitive(f.gcd(g))(1)))
            assert allowed % line == 0 and actual % allowed == 0
            count += 1
        strict += allowed > base
    print(f"Base {p}: {len(candidates)} offsets, {count} lines, {strict} larger integer divisors")


def exponential_examples():
    for length in (4, 8, 30, 64, 250, 1024):
        for p, positions in ((2, (0, 2*length, 4*length+1, length+2)),
                             (3, (0, 2*length, 4*length, length+1))):
            f = 1+X**length
            g = sum(X**j for j in positions)
            assert f.gcd(g) == 1
            assert sum(p**j for j in positions) % (1+p**length) == 0
            assert len(set(positions)) == 4
    print("Twelve exact exponential examples with trivial polynomial gcd")


if __name__ == '__main__':
    for p, (parameters, cases) in CASES.items():
        phase_controls(p, parameters, cases)
        print("Direct characteristic-zero phase orders:", dict(direct_phase_counts(p, parameters[-1], cases)))
        integer_controls(p, parameters[-1], offsets(cases))
    exponential_examples()
