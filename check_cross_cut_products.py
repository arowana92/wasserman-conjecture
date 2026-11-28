"""Enumerate endpoint-coprime cut pairs with h*L bounded. Requires python-flint."""
import argparse
from array import array
from functools import lru_cache
from math import gcd, isqrt
from flint import fmpz
from check_cross_cut_audit import common_primes


def least_prime_factors(limit):
    least = array("I", [0])*(limit+1)
    for p in range(2, isqrt(limit)+1):
        if least[p]:
            continue
        for v in range(p*p, limit+1, p):
            if not least[v]:
                least[v] = p
    return least


def root_table(limit, ell, least):
    """Use the smallest prime power and CRT to enumerate all roots modulo x."""
    @lru_cache(None)
    def local(p, power):
        if p == 2 or ell%p == 0:
            roots = [r for r in range(p) if (r*r+ell)%p == 0]
            modulus = p
            while roots and modulus < power:
                higher = modulus*p
                roots = [v for r in roots for v in range(r, higher, modulus)
                         if (v*v+ell)%higher == 0]
                modulus = higher
            return roots
        if pow(-ell, (p-1)//2, p) != 1:
            return []
        r = int(fmpz(-ell).sqrtmod(p))
        modulus = p
        inverse = pow(2*r, -1, p)
        while modulus < power:
            r += modulus*((-(r*r+ell)//modulus*inverse) % p)
            modulus *= p
        return (r, power-r)

    def roots(x):
        factors = []
        while x > 1:
            p = least[x] or x
            power = p
            x //= p
            while x%p == 0:
                power *= p
                x //= p
            choices = local(p, power)
            if not choices:
                return ()
            factors.append((power, choices))
        modulus, values = factors[0]
        for q, choices in factors[1:]:
            inv = pow(modulus, -1, q)
            values = [r+modulus*(((s-r)*inv) % q) for r in values for s in choices]
            modulus *= q
        return values

    for x in range(2, limit+1):
        for m in roots(x):
            if m:
                yield x, m
    local.cache_clear()


def scan(limit, bound):
    least = least_prime_factors(limit)
    lowest = parity = 0
    survivors = []
    for ell in range(1, bound+1):
        for x, m in root_table(limit, ell, least):
            k = m-(m*m+ell)//x
            if k < 1:
                continue
            r = gcd(m,x)
            for h in range(r, bound//ell+1, r):
                if k*m%h:
                    continue
                y, rem = divmod(m*m*(x-m)-h, ell)
                if rem or y < x or gcd(x,y) != 1 or k*x*y%h:
                    continue
                s = gcd(m,y)
                if h%gcd(r*r*(y-m),s*s*(x-m)):
                    continue
                lowest += 1
                n, g = k*x*y//h, k*m//h
                a, b = g*x, g*y
                if a & (n-a) and b & (n-b):
                    continue
                parity += 1
                if not any(common_primes(n, a, b)):
                    survivors.append((h, ell, x, m, n, a, b))
        if ell%10 == 0:
            print("Through L =", ell, "lowest layers:", lowest, "parity:", parity, flush=True)
    return lowest, parity, sorted(survivors)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("limit", type=int)
    parser.add_argument("product", type=int)
    options = parser.parse_args()
    counts = scan(options.limit, options.product)
    print("Lowest layers, parity, endpoint-coprime pairs:", *counts)
    if (options.limit, options.product) == (10**6, 50):
        assert counts == (5383328, 98, [(3, 7, 23, 4, 989, 92, 172),
                                       (10, 5, 10, 5, 46, 10, 23)])
