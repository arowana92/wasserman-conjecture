"""Cross-check the first bounded-digit obstruction using only residue sums."""

from functools import lru_cache
from math import gcd

from check_cofactor_index_depth import locally_possible
from check_prime_power_cofactors import shapes
from search_endpoint import atoms, smallest_factors


@lru_cache(maxsize=None)
def companions(u, required):
    return [(f, g, u-f-g) for f in range(1, u//3+1) for g in range(f, (u-f)//2+1)
            if all(sum(x % w == 0 for x in (f, g, u-f-g)) >= 2 for w in required)]


def independent_candidates(d, factors):
    # For a unit scalar the low residues 0,0,1 allow the unshifted pair
    # (0,0), (unit,0), or (0,unit). Otherwise C needs two zero slots.
    rows = set()
    for a in range(3, (d-1)//2+1):
        for b in range(a, d-a):
            ab = gcd(a, b)
            if ab < 3:
                continue
            required = tuple(w for r, w in factors
                             if not ((a % w == 0 and (b % w == 0 or b % r)) or
                                     (b % w == 0 and a % r)))
            for u in range(3, ab+1):
                if ab % u == 0 and gcd(u, d) == 1:
                    for f, g, h in companions(u, required):
                        rows.add((u, a, b, d-a-b, f, g, h))
    return rows


def first_empty_depth(d, shape, r, limit):
    u, a, b, c, f, g, h = shape
    # Separate row flags persist through every inspected layer on this path.
    prefixes = [(0, True, True)]
    w = 1
    for k in range(1, limit+1):
        old, w = w, r*w
        nxt = []
        inverse_u = pow(u, -1, w) if u % r else None
        for low, a_ok, c_ok in prefixes:
            for x in range(low, w, old):
                if x % r == 0:
                    continue
                n = 1+d*x
                first = a_ok and (a*x % w + b*x % w + (c*x+1) % w < w)
                if u % r == 0:
                    if gcd(n, w) != gcd(u, w) or not first:
                        continue
                    nxt.append((x, True, False))
                else:
                    v = n*inverse_u % w
                    other = c_ok and (f*v % w + g*v % w + h*v % w < w)
                    if (d % r == 0 and (first or other)) or (d % r and (n % r or other)):
                        nxt.append((x, first, other))
        if not nxt:
            return k
        prefixes = nxt
    raise ArithmeticError((d, shape, r, limit))


def verify():
    spf = smallest_factors(248)
    survivors = []
    counts = {}
    count = 0
    necessary = set()
    cover = {2, 3, 5, 7, 13, 17, 23, 47}
    for d in range(3, 248):
        factors = atoms(d, spf)
        rows = independent_candidates(d, factors)
        assert rows == set(shapes(d, factors))
        for shape in sorted(rows):
            count += 1
            primes = {r for r, _ in factors+atoms(shape[0], spf)+atoms(d+1, spf)}
            bad = {r for r in primes if not locally_possible(d, shape, 1, r)}
            if bad:
                k = first_empty_depth(d, shape, min(bad & cover), 64)
                counts[k] = counts.get(k, 0)+1
                if len(bad) == 1:
                    necessary.update(bad)
            else:
                survivors.append((d, shape))
    assert survivors == [(247, (3, 114, 117, 16, 1, 1, 1))]
    assert necessary == cover and count == 607 and sum(counts.values()) == 606
    assert max(counts) == 24
    print('Smallest possible prime set after the initial residue conditions:', sorted(cover))
    print('Candidate row shapes:', count, 'Direct exclusions by digit depth:', sorted(counts.items()))
    print('Survivors:', survivors)


if __name__ == '__main__':
    verify()
