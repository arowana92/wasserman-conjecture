"""Cross-check the exponent obstructions with pair-of-triples enumeration."""
from functools import lru_cache
from itertools import combinations_with_replacement, product
from math import gcd, isqrt

from flint import fmpz, fmpz_mod_ctx

from check_seven_perfect_powers import PRIMES


@lru_cache(None)
def prime_proof(r):
    if r < 10**6:
        assert r >= 2 and all(r % j for j in range(2, isqrt(r)+1))
        return
    remaining = r-1
    primes = []
    for p, e in fmpz(r-1).factor():
        p = int(p)
        prime_proof(p)
        power = p**e
        assert remaining % power == 0
        remaining //= power
        primes.append(p)
    assert remaining == 1
    assert any(pow(g, r-1, r) == 1 and
               all(gcd(pow(g, (r-1)//p, r)-1, r) == 1 for p in primes)
               for g in range(2, 1000))


def test_sum(ell, moduli, h):
    low = gcd(7**h-1, 432)
    small = {}
    large = []
    for r in moduli:
        assert pow(7, h, r) == 1 and (r-1) % ell == 0
        if r < 150000:
            small[r] = {(6*pow(y, ell, r)+6) % r for y in range(r)}
        else:
            ring = fmpz_mod_ctx(r)
            k = (r-1)//ell
            large.append((r, ring, k, ring(6)**k))

    def accepts(n):
        if (n-72) % low or any(n % r not in values for r, values in small.items()):
            return False
        # For n = 6*y^ell + 6, a nonzero (n-6)^((r-1)/ell)
        # agrees with 6^((r-1)/ell) by Fermat.
        for r, ring, k, target in large:
            x = (n-6) % r
            if x and ring(x)**k != target:
                return False
        return True

    return accepts


def seeds(ell, accept):
    w = [7**i for i in range(ell)]
    triples = [(c, sum(w[i] for i in c))
               for c in combinations_with_replacement(range(ell), 3)]
    # Each sorted sextuple splits uniquely after its third term.
    heads = [[(xs, s) for xs, s in triples if xs[-1] <= i] for i in range(ell)]
    result = set()
    for tail, upper in triples:
        for head, lower in heads[tail[0]]:
            if accept(lower+upper):
                result.add(head+tail)
    return result


def check(ell):
    pending = set(PRIMES[ell])
    for r in pending:
        prime_proof(r)
    live = set()
    counts = []
    old_h = 1
    for factor in (1, 2, 6, 12, 36):
        h = factor*ell
        active = {r for r in pending if pow(7, h, r) == 1}
        pending -= active
        accept = test_sum(ell, sorted(active), h)
        if factor == 1:
            live = seeds(ell, accept)
        else:
            following = set()
            powers = [7**i for i in range(h)]
            for previous in live:
                choices = [range(i, h, old_h) for i in previous]
                for xs in product(*choices):
                    if accept(sum(powers[i] for i in xs)):
                        following.add(tuple(sorted(xs)))
            live = following
        counts.append(len(live))
        if not live:
            print(ell, counts, flush=True)
            return
        old_h = h
    raise ValueError((ell, "surviving tuples", len(live)))


if __name__ == "__main__":
    from sys import argv
    for exponent in map(int, argv[1:]) if argv[1:] else PRIMES:
        check(exponent)
