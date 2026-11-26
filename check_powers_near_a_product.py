"""Check the 122,070 multipliers for a prime power close to a product."""
from itertools import combinations

from check_quartic_triangle import binomial_order
from check_scaled_row_blocks import certify_prime
from check_signed_power_pairs import pair_test
from search_endpoint import atoms, smallest_factors


def small_examples():
    table = smallest_factors(32400)
    powers = {v: factors[0][0] for v in range(2, 32401)
              if len(factors := atoms(v, table)) == 1}
    checked = locked = carries = 0
    for p, q in combinations(sorted(v for v in powers if v <= 180), 2):
        if powers[p] == powers[q]:
            continue
        u = q-p
        if u > 12:
            continue
        for k in range(1, min(u*u, p, 9)):
            r = p*q-k
            if r not in powers or powers[r] in (powers[p], powers[q]):
                continue
            for d in range(1, 6):
                if d*k >= p:
                    continue
                n = d*p*q*r
                for w in (p, q):
                    for cut in range(w*r, n, w*r):
                        checked += 1
                        if binomial_order(n, cut, powers[r]):
                            continue
                        assert cut % (p*q*r) == 0
                        locked += 1
                        if 4*(d*max(k*u, u*u, 2*u)+1) < w:
                            assert binomial_order(n, cut, powers[w]) > 0
                            carries += 1
    assert 0 < carries <= locked < checked
    print(f"{checked} binomial cut tests, {locked} forced multiples, {carries} digit exclusions.")


def example():
    p, q, r = 100000000003, 100000000019, 10000000002200000000023
    u, k = q-p, p*q-r
    for w in (p, q, r):
        certify_prime(w)
    f = p*q*r
    limit = p//(200*u**3)
    assert (u, k, limit) == (16, 34, 122070) and k < u*u
    for d in range(1, limit+1):
        assert pair_test(d*f, p, q)
        assert d*k < p and d*k*u < p and d*(2*u) < q
    print(f"All three primes certified and all multipliers through {limit} checked.")


if __name__ == "__main__":
    small_examples()
    example()
