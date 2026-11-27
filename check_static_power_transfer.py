"""Verify the cross-prime congruence and count reference cuts for four prime powers."""

from argparse import ArgumentParser
from math import gcd

from check_two_power_seed import free, grouped, split_costs
from prime_certificates import primes_through


def prime_powers(limit):
    result = {}
    for p in primes_through(limit):
        w = p
        while w <= limit:
            result[w] = p
            w *= p
    return result


def allowed(a, b, p):
    if b % p == 0 or a <= 3:
        return False
    if a == 4:
        return b % 3 == 2 and ((b+1) & (3*b)) == 0
    n = a*b+1
    for c in range(2, a):
        if free((0, c*b, n-c*b), p):
            return True
    for w in range(2, a):
        if n % w:
            continue
        h = n//w
        if any(free((0, f*h, (w-f)*h), p) for f in range(1, w//2+1)):
            return True
    return False


def check_gluing(limit):
    powers = prime_powers(limit)
    candidates = equalities = comparisons = 0
    for a, p in powers.items():
        for b in powers:
            d = a*b
            if d > limit or b % p == 0:
                continue
            n = d+1
            if not allowed(a, b, p):
                assert split_costs(d)[1] >= d and next(grouped(d), None) is None, (a, b)
                comparisons += 1
            for x in range(a, d//2, a):
                for y in range(b, (d-1)//2-x+1, b):
                    if gcd(x, y) != 1:
                        continue
                    candidates += 1
                    at_a = pow(y, -1, a)
                    at_b = pow(x, -1, b)
                    if max(at_a, at_b)*(x+y) <= n:
                        assert at_a == at_b and at_a*(x+y) == n, (a, b, x, y)
                        equalities += 1
    print(f"Through cofactor {limit}: {candidates} primitive pairs, {equalities} "
          f"exact equalities, {comparisons} exclusions cross-checked by the general test.")


def affine_prefixes():
    checked = delayed_total = 0
    for p in (2, 3, 5):
        for a in range(1, 6):
            for b in range(1, 6):
                valuations = []
                for c in (a, b, a+b):
                    v = 0
                    while c % p == 0:
                        c //= p
                        v += 1
                    valuations.append(v)
                base = min(valuations)
                delay = max(valuations)-base
                for m in (3, 4):
                    bound = 1
                    if delay == 0:
                        bound = ((p+1)//2)**m
                    else:
                        for start in range(min(delay, m)):
                            chains = [1]*p
                            for _ in range(start+delay, m, delay):
                                chains = [sum(chains[:p-j]) for j in range(p)]
                            bound *= sum(chains)
                    for i in (0, 1):
                        for j in (0, 1):
                            accepted = 0
                            for z in range(p**m):
                                x, y = a*z+i, b*z+j
                                for _ in range(m+base):
                                    x, digit = divmod(x, p)
                                    y, other = divmod(y, p)
                                    if digit+other >= p:
                                        break
                                else:
                                    accepted += 1
                            assert accepted <= bound, (p, a, b, m, i, j, accepted, bound)
                            checked += 1
                            delayed_total += valuations[2] > max(valuations[:2])
    print(f"{checked} exact affine prefix counts, {delayed_total} with the total delayed.")


def census(limit):
    powers = prime_powers(limit)
    for a in (2, 3, 4, 5, 7, 8, 9, 11, 13, 16):
        p = powers[a]
        eligible = sum(not allowed(a, b, p) for b in powers)
        print(f"a={a}: {eligible} of {len(powers)} prime powers b through {limit} certified.")


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--cross", type=int, default=500)
    parser.add_argument("--limit", type=int, default=100000)
    args = parser.parse_args()
    if args.cross < 1 or args.limit < 16:
        parser.error("need a positive cross-check range and limit >= 16")
    check_gluing(args.cross)
    affine_prefixes()
    census(args.limit)
