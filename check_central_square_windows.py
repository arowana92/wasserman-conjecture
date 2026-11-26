"""Find every unreflected cut pair in a central window that passes the endpoint tests."""
from argparse import ArgumentParser
from collections import Counter, defaultdict
from itertools import combinations_with_replacement
from math import comb, gcd, isqrt

from check_central_cuts import (binomial_carry, endpoint_primes, legendre,
                                lower_powers_pass, prime_list, quadratic_candidates,
                                refinement_masks)


def eligible(n, x, y, c, power):
    if not 0 <= y < x <= n-2 or x*x > c*n or power-n > x:
        return None
    if (n-x) % 2 or (n-y) % 2:
        return None
    a, b = (n-x)//2, (n-y)//2
    if (a & (n-a)) and (b & (n-b)):
        return None
    if lower_powers_pass(n, a, b):
        return n, a, b
    return None


def by_square(c, bits, difference=False):
    """Sweep intervals of squares for x+y or x-y. In the latter sweep,
    (x-y)^2=(t-1)^2 uses the earlier fixed-gap reconstruction separately."""
    if difference:
        # Above this point any nonzero coefficient of M in (x-y)^2 is positive.
        low = (c+1)**2+2
        yield from direct(c, min(low-1, (1 << bits)-1))
        for n, a, b in quadratic_candidates((c-1)//2, bits):
            if n >= low and (n-2*a)**2 <= c*n and lower_powers_pass(n, a, b):
                yield n, a, b
    else:
        low = 3
    sign = -1 if difference else 1
    for k in range(1, 4*(c+1)):
        targets = defaultdict(list)
        for t in range(c+1):
            targets[(t+sign)**2 % k].append((t, (t+sign)**2))
        for m in range(max(2, low.bit_length()), bits+1):
            power = 1 << m
            bot = max(low, power//2, power-isqrt(c*(power-1)))
            if bot >= power:
                continue
            first = isqrt(k*(bot-1))
            last = isqrt(k*(power-2)+(c+1)**2)
            for v in range(first, last+1):
                square = v*v
                for t, constant in targets.get(square % k, ()):
                    n = (square-constant)//k+1
                    if not bot <= n < power:
                        continue
                    other_square = square-sign*4*t*n
                    if other_square < 0 or other_square % 16 not in (0, 1, 4, 9):
                        continue
                    if other_square % 63 not in SQUARES_63:
                        continue
                    w = isqrt(other_square)
                    if w*w != other_square or (v+w) % 2:
                        continue
                    x, y = (v+w)//2, sign*((v-w)//2)
                    match = eligible(n, x, y, c, power)
                    if match:
                        yield match


SQUARES_63 = {i*i % 63 for i in range(63)}


def direct(c, last):
    for n in range(3, last+1):
        power = 1 << n.bit_length()
        radius = isqrt(c*n)
        if power-n > radius:
            continue
        indices = range(max(1, (n-radius+1)//2), n//2+1)
        # A cut pair passing the binary test contains at least one odd cut.
        odd = [a for a in indices if (a & (n-a)) == 0]
        found = set()
        for a in odd:
            for b in indices:
                i, j = sorted((a, b))
                if i != j and lower_powers_pass(n, i, j):
                    found.add((n, i, j))
        yield from found


def verify(c, bits, cross_bits):
    found = set(by_square(c, bits))
    control = {case for case in found if case[0].bit_length() <= cross_bits}
    assert control == set(by_square(c, cross_bits, difference=True))
    small = min((1 << bits)-1, 12000)
    assert {case for case in found if case[0] <= small} == set(direct(c, small))
    primes = prime_list(5000)
    witnesses = Counter()
    exceptions = []
    for n, a, b in sorted(found):
        for p in primes:
            if n % p <= 1 and binomial_carry(n, a, p) and binomial_carry(n, b, p):
                assert all(legendre(n, p) > legendre(i, p)+legendre(n-i, p) for i in (a,b))
                witnesses[p] += 1
                break
        else:
            if n > 1_000_000:
                raise ArithmeticError(f"Unresolved candidate of bit length {n.bit_length()}.")
            common = gcd(n*(n-1), gcd(comb(n, a), comb(n, b)))
            if common > 1:
                witnesses[common] += 1
                continue
            ps = endpoint_primes(n)
            ma, mb = (refinement_masks(n, i, ps) for i in (a, b))
            assert all(left & right for left in ma for right in mb)
            exceptions.append((n, a, b))
    print(f"Squared offset limit {c}N, weights below 2**{bits}.")
    print(f"{len(found)} candidates, {len(control)} in the second sweep through {cross_bits} bits.")
    print(f"Endpoint witnesses: {dict(witnesses)}.")
    print(f"All trinomial refinements covered at the binomial exceptions: {exceptions}.")


def algebra_controls():
    for n in range(3, 151):
        bins = [comb(n, a) for a in range(n//2+1)]
        for a, b in combinations_with_replacement(range(1, n//2+1), 2):
            x, y = n-2*a, n-2*b
            common = gcd(n*(n-1), gcd(bins[a], bins[b]))
            if x*x <= 5*n:
                assert common > 1
            if common != 1:
                continue
            t, rem = divmod(x*y, n)
            assert rem == 0 and t != 1
            k, rem = divmod((x+y)**2-(t+1)**2, n-1)
            assert not rem and (x-y)**2 == (k-4*t)*(n-1)+(t-1)**2
    for y in range(2, 1000):
        m = 1
        while m < y:
            assert y*(y-m) & (y*(y+m))
            m *= 2


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--limit", type=int, default=64)
    parser.add_argument("--bits", type=int, default=4096)
    parser.add_argument("--cross-bits", type=int, default=256)
    args = parser.parse_args()
    if args.limit < 1 or not 2 <= args.cross_bits <= args.bits:
        parser.error("Require a positive limit and 2 <= cross-bits <= bits.")
    algebra_controls()
    verify(args.limit, args.bits, args.cross_bits)
