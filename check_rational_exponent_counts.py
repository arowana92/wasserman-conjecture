"""Exact controls for the rational-function and power-tower sampling lemmas."""
from fractions import Fraction
from itertools import combinations
from math import inf, lcm


def valuation(x, p):
    x = Fraction(x)
    if not x:
        return inf
    n, d = abs(x.numerator), x.denominator
    value = 0
    while n % p == 0:
        n //= p
        value += 1
    while d % p == 0:
        d //= p
        value -= 1
    return value


def rational_balls():
    total = 0
    # R = c*(X-r)^a/(X-s)^b. Its distinct finite derivative zeros
    # and poles can be found explicitly even for large multiplicities.
    for r, s, a, b, c in [(Fraction(1,2), Fraction(-3,5), 27, 5, Fraction(7,16)),
                           (Fraction(-4), Fraction(3), 13, 13, Fraction(1,27)),
                           (Fraction(3,5), Fraction(-4), 1, 8, Fraction(5,49))]:
        centers = {s}
        if a > 1:
            centers.add(r)
        if a != b:
            centers.add((a*s-b*r)/(a-b))
        values, derivatives = {}, {}
        for x in range(-48, 65):
            if x == s:
                continue
            values[x] = c*(x-r)**a/(x-s)**b
            derivatives[x] = c*(x-r)**(a-1)*((a-b)*x+b*r-a*s)/(x-s)**(b+1)
        for p in (2,3,5,7,11):
            heights = {x: max(valuation(x-u,p) for u in centers) for x in values}
            for x, y in combinations(values, 2):
                depth = valuation(y-x, p)
                if heights[x] < depth-1:
                    assert valuation(values[y]-values[x], p) == depth + valuation(derivatives[x], p)
                    total += 1
            for n in range(1, 5):
                near = {x % (p**n) for x in values if heights[x] >= n-1}
                assert len(near) <= p*len(centers)
    # Nonsplit quadratic derivative X^2-p. Both roots have valuation 1/2.
    for p in (2,3,5,7):
        for x, y in combinations(range(120), 2):
            n = valuation(y-x,p)
            nearest = Fraction(1,2) if x % p == 0 else 0
            if nearest < n-1:
                delta = Fraction(y**3-x**3,3)-p*(y-x)
                assert valuation(delta,p) == valuation(x*x-p,p)+n
                total += 1
    return total


def tower_checks():
    comparisons = 0
    for p in (2,3,5,7,11):
        phase = lcm(*range(1,p)) if p > 2 else 2
        cutoff = max(1, phase.bit_length())
        for a in (3,4,5,7,10):
            for b in (3,4,7,9):
                if a % p == 0 or b % p == 0:
                    continue
                def gain(base):
                    k = 2 if p == 2 else p-1
                    return valuation(base**k-1,p)-valuation(k,p)
                for j in (1,2,p,2*p,p*p):
                    x = cutoff+1
                    y = x+phase*j
                    predicted = valuation(y-x,p)+gain(a)+gain(b)
                    modulus = p**(predicted+1)
                    # Reduce each exponent using Euler only on the outer unit.
                    phi = (p-1)*(modulus//p)
                    left = pow(b,pow(a,x,phi),modulus)
                    right = pow(b,pow(a,y,phi),modulus)
                    assert valuation(left-right,p) == predicted
                    assert pow(a,x,phase) == pow(a,y,phase)
                    comparisons += 1
    return comparisons


def prefix_control():
    # 9^j ranges over every unit congruent to 1 modulo 8.
    counts = [1,2,3]
    for n in range(3,15):
        counts.append(counts[n-1]+counts[n-3])
    for k, expected in enumerate(counts):
        modulus = 2**(k+3)
        v, good = 1, 0
        for _ in range(2**k):
            a, b, c = v, 2*v, 4*v
            if ((a&b | a&c | b&c) % modulus) == 0:
                good += 1
            v = 9*v % modulus
        assert v == 1 and good == expected
    upper = Fraction(1466,1000)
    assert upper**3 > upper**2+1 and upper**9 < 32
    return sum(counts)


if __name__ == "__main__":
    print(rational_balls(), "local divided-difference identities")
    print(tower_checks(), "two-layer valuation comparisons")
    print(prefix_control(), "accepted low prefixes across fifteen orbit lengths")
