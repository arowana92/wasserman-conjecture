"""Verify a shared minimum certificate using exact binomial valuations."""
from itertools import combinations, product
from math import comb, gcd, prod

from check_scaled_row_blocks import certify_prime


def digits(value, base):
    out = []
    while value:
        value, digit = divmod(value, base)
        out.append(digit)
    return out


def binomial_order(n, cut, prime):
    total, first, second = n, cut, n-cut
    order = 0
    while total:
        total //= prime
        first //= prime
        second //= prime
        order += total-first-second
    return order


def common_cuts(n, w, z, p, q):
    """All proper cuts avoiding p and q for w and z powers of the bases."""
    assert w < z and n % (w*z) == 0 and gcd(w, z) == 1
    bounds = digits(n//w, w)
    place = [w**i for i in range(1, len(bounds))]
    tested = 0
    avoid = set()
    for high in product(*(range(t+1) for t in bounds[1:])):
        upper = sum(x*y for x, y in zip(high, place))
        low = (-upper) % z
        if low > bounds[0]:
            continue
        cut = w*(upper+low)
        if 0 < cut < n:
            tested += 1
            if binomial_order(n, cut, p) == 0 and binomial_order(n, cut, q) == 0:
                avoid.add(cut)
    return tested, avoid


def cross_check():
    count = allowed = excluded = 0
    powers = ((2, 2), (3, 3), (4, 2), (5, 5), (7, 7), (8, 2), (9, 3))
    for (w, p), (z, q) in combinations(powers, 2):
        if p == q:
            continue
        for m in range(1, 15):
            n = w*z*m
            # Every coefficient is formed as an exact integer in this control.
            direct = {k for k in range(1, n) if gcd(comb(n, k), p*q) == 1}
            _, found = common_cuts(n, w, z, p, q)
            assert direct == found, (w, z, n, direct ^ found)
            count += 1
            allowed += len(found)
            excluded += not found
    assert allowed and excluded
    print(f"{count} exact small comparisons, {allowed} avoiding cuts, {excluded} empty sets.")


def certificate():
    p = 1389209
    a, b, c = 2, 8, 24
    s, t = b-a, c-b
    weights = [p+h for h in (0, a, b, c)]
    n = prod(weights)
    r = n-2
    assert p >= 100*c**3 and s*(s+t) > a*(2*s+t) and t > b+s
    assert n == 3724611517745287445819539 and n+1-r <= p and 2*r > n+1
    factors = [(p0, 1) for p0 in weights]
    successor = [(2, 2), (5, 1), (42683, 1), (4363108869743560019, 1)]
    for q in {r} | {q for q, _ in factors+successor}:
        certify_prime(q)
    for weight, decomposition in ((n, factors), (n+1, successor)):
        assert prod(q**v for q, v in decomposition) == weight
        assert all(sum(digits(weight, q)) >= 3 for q, _ in decomposition)
    # For a prime not dividing a weight, digit sum <= 2 would force
    # weight-1 to be a pure prime power. Both predecessors have two bases.
    assert (n-1) % 6 == 0 and n % (p*(p+a)) == 0
    for w, z in combinations(weights[:3], 2):
        tested, found = common_cuts(n, w, z, w, z)
        assert not found
        print(f"{w}, {z}: {tested} proper cuts passed the first block test, zero avoided both.")
    print(f"Minimum certificate at {n} and {n+1}: {weights[:3] + [r]}")


if __name__ == "__main__":
    cross_check()
    certificate()
