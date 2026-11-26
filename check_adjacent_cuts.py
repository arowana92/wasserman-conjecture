#!/usr/bin/env python3
"""Adjacent binomial cuts and the prime factors of N(N-1)."""
import argparse
from math import comb, gcd, isqrt, lcm

from check_central_cuts import binomial_carry, endpoint_primes


def check_binary_base():
    allowed = []
    for d in range(0, 16, 2):
        n = (d*d - 1) * pow(3, -1, 16) % 16
        a = (n - d - 1) // 2
        if any(((c % 8) & ((n-c) % 8)) == 0 for c in (a, a+1)):
            allowed.append(d)
    assert allowed == [2, 4, 12, 14]
    for r in (2, 4, -4, -2):
        n = (r*r - 1) // 3
        a = (n-r-1)//2
        for c in (a, a+1):
            assert any(0 <= side <= 5 for side in (c, n-c))


def check_edge_parameters(g, ell, i, j, c):
    n, a = 1+g*i, g*c
    k = c-g*j
    edge = k-j*ell
    assert 0 < k < 4*j*ell and 0 < edge*g < 12*j*ell**2
    w = edge*g-ell*k
    assert j*w*w == ell*k**3-edge**2
    z = 2*w-2*ell**2*j-3*ell*edge
    assert 0 < 4*j*z < 3*edge**2
    t = 3*edge**2-4*j*z
    assert t*(2*ell**2*j+3*ell*edge) == ell*edge**3+2*j*z*z+8*edge**2
    v = j*ell**2
    assert 108*a**5 > z*(n-1)**4
    residual = (3*g*g-14*g*ell+14*ell**2)*edge - v*(6*g-13*ell)
    if g >= 16*ell:
        assert abs(residual)*g**3 < 2000*ell**4*v
        assert 0 < t*g**3 < 12100*v*v*ell
        if residual:
            assert 470**4*n**4*a**21 > (n-1)**21
    return residual


def zero_residuals():
    cases = set()
    for ell in (1, 2, 3, 6, 9, 18):
        multiple = 354294*ell**3
        for d in range(1, isqrt(multiple)+1):
            if multiple % d:
                continue
            for divisor in {d, multiple//d}:
                g, rem = divmod(divisor-ell, 12)
                if rem or g <= 4*ell or gcd(g, ell) != 1:
                    continue
                j, rem = divmod((6*g-13*ell)*(3*g*g-14*g*ell+14*ell**2),
                                ell**4*divisor)
                if not rem and j > 0:
                    cases.add((g, ell, j))
    assert cases == {(20, 3, 2), (182, 3, 546)}
    assert (1+182**2*546) % 3 and 20 < 16*3


def enumerate_pairs(limit):
    tested = 0
    low = set()
    for g in range(5, isqrt(limit)//2 + 1):
        for ell in range(1, (g-1)//4 + 1):
            if gcd(g, ell) != 1:
                continue
            first = (-pow(g*g, -1, ell)) % ell if ell > 1 else 1
            last = ((limit-1)*ell-g) // (g**3)
            for j in range(first, last+1, ell):
                tested += 1
                i = (1+g*g*j)//ell
                h2 = i*(i-4*g*j)-4*j
                if h2 < 0:
                    continue
                h = isqrt(h2)
                if h*h != h2 or (i-h) % 2:
                    continue
                c, f = (i-h)//2, (i+h)//2
                if c <= 0 or gcd(c, f) != 1:
                    continue
                n, a = 1+g*i, g*c
                if 2*(a+1) > n:
                    continue
                assert gcd(a, n-1) == g and c*f == j*n
                d = n-2*a-1
                delta = g-4*ell
                e, rem = divmod(1+delta*delta*j, ell)
                assert rem == 0 and d*d == 3*n+1+n*i*delta
                assert 9*d**6 > 4*n**4*(n-1)*delta
                z = 4*delta*j + e
                assert delta*e*h*h == delta*z**3+(z-e)*(4*z-e)
                check_edge_parameters(g, ell, i, j, c)
                low.add((n, a))
    return tested, low


def direct_low(limit):
    """Scan indices independently. Check each full weight or predecessor power."""
    cases = set()
    for n in range(3, limit+1):
        powers = [[], []]
        for p in endpoint_primes(n):
            kind = 0 if n % p == 0 else 1
            w = p
            while (n-kind) % (w*p) == 0:
                w *= p
            powers[kind].append(w)
        for a in range(1, n//2):
            if 4*(a+1)*(n-a) >= n*(n-1):
                continue
            if all(a % w in (0, w-1) for w in powers[0]) and all(
                    a % w in (0, 1, w-1) for w in powers[1]):
                cases.add((n, a))
    return cases


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--limit", type=int, default=10**9)
    parser.add_argument("--direct", type=int, default=4096)
    args = parser.parse_args()
    if args.limit < 3 or not 3 <= args.direct <= args.limit:
        parser.error("require 3 <= direct <= limit")
    check_binary_base()
    zero_residuals()
    tested, algebraic = enumerate_pairs(args.limit)
    low = {(n, a) for n, a in algebraic if lcm(a-1, a, a+1) % (n-1) == 0}
    assert {pair for pair in low if pair[0] <= args.direct} == direct_low(args.direct)
    exact_count = 0
    for n in range(3, min(args.direct, 512)+1):
        for a in range(1, (n-1)//2 + 1):
            exact_count += 1
            assert gcd(n*(n-1), gcd(comb(n, a), comb(n, a+1))) > 1
    live = [(n, a) for n, a in sorted(algebraic)
            if not (binomial_carry(n, a, 2) and binomial_carry(n, a+1, 2))]
    unresolved = []
    for n, a in live:
        if not any(binomial_carry(n, a, p) and binomial_carry(n, a+1, p)
                   for p in endpoint_primes(n)):
            unresolved.append((n, a))
    print(f"weights through {args.limit}, parameter triples {tested}, algebraic pairs {len(algebraic)}, low-power pairs {len(low)}")
    print(f"exact adjacent coefficient checks {exact_count}, binary survivors {len(live)}")
    print("endpoint gcd one", unresolved)


if __name__ == "__main__":
    main()
