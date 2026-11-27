"""Check polynomial lifts, delayed odd-prime digits, and powers as input parameters."""

from fractions import Fraction
from itertools import combinations, product
from math import comb

from check_endpoint_pair_bound import prime_powers


def value(f, x):
    total = 0
    for coefficient in reversed(f):
        total = total*x+coefficient
    return total


def valuation(n, p):
    assert n
    e = 0
    while n % p == 0:
        n //= p
        e += 1
    return e


def signed_controls():
    # For each trial rational root, test the homogeneous sum against exact
    # rational evaluation and integer synthetic division.
    checked = roots = negative = 0
    for a, b, c in product(range(1, 5), range(-4, 5), range(1, 6)):
        f = (1, a, b, c)
        for x in range(1, 8):
            for k in range(1, 10):
                norm = sum(v*(-x)**i*k**(3-i) for i, v in enumerate(f))
                assert norm == k**3 * value(f, Fraction(-x, k))
                checked += 1
                if norm:
                    continue
                roots += 1
                t, rem = divmod(k, x)
                assert not rem
                g = (1, a-t, b-t*(a-t))
                assert t*g[-1] == c
                if min(g) >= 0:
                    continue
                negative += 1
                # Test the first negative column against the actual two
                # digit streams in an ordinary prime radix.
                place = next(i for i, v in enumerate(g) if v < 0)
                P = 257
                for i in range(1, a):
                    j = a-i
                    numbers = (i*value(g, P), j*value(g, P))
                    assert sum(n//P**place % P for n in numbers) > P
    d, P, R = 119, 823663, 34206080687
    f = (1, 2, -1, 6)
    h = f[1]
    bound = max(abs(sum(v*(-x)**i*k**(3-i) for i,v in enumerate(f)))
                for x in range(1, (d+1)//2) for k in range(1, (d*h+3)//4))
    assert bound == 1643032 and 2*P > bound
    assert all(list(prime_powers(n)) == [(n,n)] for n in (P, R))
    assert value(f, P) == 1+d*P*R == (1+3*P)*(1-P+2*P*P)
    assert roots and negative
    print(f"Homogeneous evaluations: {checked}, rational roots: {roots}, "
          f"signed quotient controls: {negative}. Example weight: {value(f,P)}.")


def digit_blocks():
    prefixes = tested = passing = 0
    for r in (3, 5, 7):
        length = 3
        for a, b in combinations(range(1, 11), 2):
            v = min(valuation(a,r), valuation(b,r))
            delay = max(valuation(n,r) for n in (a,b,a+b))-v
            for c, d in ((0,0), (1,0), (12,9), (57,23)):
                kept = 0
                for w in range(r**length):
                    first, second = a*w+c, b*w+d
                    # Legendre's criterion, only at the low prime-power layers.
                    if all(first % m + second % m < m
                           for m in (r**j for j in range(1, v+length+1))):
                        kept += 1
                initial = Fraction(4,3)**(2*delay) if r == 3 else 3**delay
                assert kept <= initial*Fraction(5*r,6)**length
                prefixes += 1
                tested += r**length
                passing += kept
        # Cross-check the floor test at all digit positions against a
        # binomial coefficient on smaller inputs.
        for first in range(25):
            for second in range(25):
                m = r
                carry = False
                while m <= first+second:
                    carry |= first % m + second % m >= m
                    m *= r
                assert carry == (comb(first+second, first) % r == 0)
    print(f"Odd-prime blocks: {prefixes} settings, {tested} words, {passing} passing.")


def local_distances():
    compared = good = 0
    for r, p in ((2,3), (2,5), (2,17), (3,2), (5,2)):
        phase = 2 if r == 2 else r-1
        for step in (1, 3, r*r):
            inputs = [p**(3+step*phase*j) for j in range(15)]
            c = valuation(inputs[1]-inputs[0], r)
            for j,k in combinations(range(len(inputs)), 2):
                x,y = inputs[j],inputs[k]
                dist = valuation(y-x, r)
                assert dist == c+valuation(k-j, r)
                compared += 1
                # S(T)=(T-b)^m, with increasing critical multiplicity.
                for b, m in ((2,2), (1,3), (4,5), (7,6)):
                    if dist <= 1+valuation(x-b, r):
                        continue
                    derivative = m*(x-b)**(m-1)
                    assert valuation((x-b)**m-(y-b)**m, r) == valuation(derivative,r)+dist
                    good += 1
    print(f"Exponential input distances: {compared}, polynomial comparisons on good balls: {good}.")


if __name__ == "__main__":
    signed_controls()
    digit_blocks()
    local_distances()
