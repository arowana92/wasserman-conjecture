"""One perfect power at nested prime powers, using two independent residue tests."""

from functools import lru_cache
from math import gcd

from check_static_power_transfer import prime_powers


@lru_cache(None)
def images(e, modulus):
    return {pow(z, e, modulus) for z in range(modulus) if gcd(z, modulus) == 1}


@lru_cache(None)
def minimal_scales(limit):
    # A scale at least three contains an odd prime or a factor four.
    scales = {p for p in prime_powers(limit).values() if p != 2}
    if limit >= 4:
        scales.add(4)
    return tuple(sorted(scales))


def possible(b, h, exponent, side):
    D = h*b
    for u in minimal_scales(D-2):
        if gcd(D, u) != 1:
            continue
        for k in range(1, D-u):
            if gcd(k, b*u) != 1:
                continue
            if side == 'P':
                if ((-pow(k, -1, b)) % b in images(exponent, b)
                        and (-pow(k, -1, u)) % u in images(exponent, u)):
                    return u, k, u
            else:
                if (k*pow(D, -1, u)) % u not in images(exponent, u):
                    continue
                for x in range(u, D-k, u):
                    if gcd(x, b) == 1 and (-k*pow(x, -1, b)) % b in images(exponent, b):
                        return x, k, u
    return None


def direct(b, h, e, side):
    """Enumerate arbitrary scale divisors and actual powers in combined residue rings."""
    D = h*b
    for x in range(3, D-1):
        if gcd(x, b) != 1:
            continue
        for u in range(3, x+1):
            if x % u or gcd(D, u) != 1:
                continue
            for z in images(e, b*u):
                if side == 'P':
                    k = -pow(z, -1, b*u) % (b*u)
                    if x+k < D:
                        return True
                else:
                    k = (-x*z) % b
                    while x+k < D:
                        if (D*z-k) % u == 0 and gcd(k, u) == 1:
                            return True
                        k += b
    return False


def check(limit=500):
    for h in (1, 2):
        for e in (2, 3, 4, 6):
            for side in ('P', 'Q'):
                certified = []
                for b, p in sorted(prime_powers(limit).items()):
                    if h == 2 and p == 2:
                        continue
                    witness = possible(b, h, e, side)
                    if witness is None:
                        certified.append(b)
                    if b <= 40 or witness is None:
                        assert bool(witness) == direct(b, h, e, side), (b, h, e, side)
                print(h, e, side, certified, flush=True)


def ten_classes():
    target = {3: {(3, 14, 13), (3, 7, 11), (3, 11, 7), (6, 14, 4), (6, 7, 8)},
              7: {(7, 34, 12), (7, 17, 24)}}
    for u, triples in target.items():
        modulus = 5*u
        found = set()
        for P in range(modulus):
            for Q in range(modulus):
                if (10*P*Q+1) % u:
                    continue
                for x in range(u, 10, u):
                    if x*P*Q % 5 == 1 and any((k*P+1) % modulus == 0
                                               for k in range(1, 10-x)):
                        found.add((x, P, Q))
        assert found == triples
    units = [r for r in range(2205) if gcd(r, 2205) == 1]
    remaining = 0
    masks = {u: {(p, q) for _, p, q in triples} for u, triples in target.items()}
    for P in units:
        for Q in units:
            n = 1+10*P*Q
            if any((P % (5*u), Q % (5*u)) in pairs and n % (u*u) != 0
                   for u, pairs in masks.items()):
                remaining += 1
    assert remaining == 55944 and len(units)**2 == 1016064
    print(f"At ten: {remaining} surviving unit residue pairs of {len(units)**2}")


if __name__ == '__main__':
    check()
    ten_classes()
