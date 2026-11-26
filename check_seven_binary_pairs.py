"""Finite binary gcd exclusion when the units-place digit belongs to the triple."""

from collections import Counter
from math import gcd

from flint import fmpz


def circle_window(n, g):
    size = (n + 2*(g+1))//3
    left = max(0, n-2*size)
    right = min(size, n-1)
    return {r % n for t in range(left, right+1) for r in (t, -t)}, size


def candidates(n, g):
    m = n//g
    window, size = circle_window(n, g)
    for k in range(1, (3*m)//2+1, 2):
        h = g*k
        if gcd(k, m) != 1 or h % n not in window:
            continue
        for i in window:
            j = (i+h) % n
            small, large = sorted((i, j))
            if max(small, large-small, n-large) <= size:
                yield h, i
                yield h, i+n


def scan():
    counts = Counter()
    for n in range(1, 2701):
        f = (fmpz(1) << n)+1
        for g in range(1, min(n//3, 28)+1):
            if n % g or (n//g) % 2 == 0:
                continue
            for h, r in candidates(n, g):
                d = (1 << g)+1
                common = f.gcd(1 + (fmpz(1) << r) * ((fmpz(1) << h)+1))
                assert common**3 * d*d < 27 << (2*max(n, h)), (n, h, r)
                counts['residues'] += 1
                counts['nontrivial gcds'] += common > 1
    print('Nontrivial binomial factor through n=2700:', dict(counts))


def direct_control():
    counts = Counter()
    for n in range(1, 86):
        f = (1 << n)+1
        powers = [pow(2, r, f) for r in range(2*n)]
        for h in range(1, 3*n//2+1):
            d = gcd(f, (1 << h)+1)
            g = gcd(n, h)
            tested = (set(candidates(n, g)) if d > 1 and n//g > 1 else set())
            for r, q in enumerate(powers):
                u = gcd(f, 1+q*((1 << h)+1))
                if u**3*d*d >= 27 << (2*max(n, h)):
                    assert d == 1 and u == f
                    assert ((h, r) == (n+1, 0) or (h, r) == (n-1, n+1))
                    counts['aligned survivors'] += 1
                if d > 1 and n//g > 1:
                    a, b = sorted((r % n, (r+h) % n))
                    assert ((h, r) in tested) == (3*max(a, b-a, n-b) <= n+2*(g+1))
                counts['unfiltered residues'] += 1
    assert min(counts.values()) > 0
    print('Small unfiltered control:', dict(counts))


def order_cross_check():
    counts = Counter()
    for n in range(1, 2701):
        f = (fmpz(1) << n)+1
        for h in range(1, (3*n)//2+1):
            g = gcd(n, h)
            if g > 28 or g == n or n//g % 2 == 0 or h//g % 2 == 0:
                continue
            d = (1 << g)+1
            cofactor = f//d
            z = (fmpz(1) << h)+1
            cap = cofactor.gcd(pow(z, 2*n, cofactor)-1)
            counts['gap pairs'] += 1
            limit = 27 << (2*max(n, h))
            if cap**3*d*d < limit:
                continue
            counts['full period checks'] += 1
            residue = z % cofactor
            for _ in range(2*n):
                u = cofactor.gcd(residue+1)
                assert u**3*d*d < limit
                residue = (2*residue) % cofactor
    print('Independent order bound:', dict(counts))


if __name__ == '__main__':
    direct_control()
    scan()
    order_cross_check()
