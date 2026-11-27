"""Finite controls for the counting theorem with two variable prime powers."""
from itertools import combinations, product
from math import comb, gcd, isqrt

from check_endpoint_pair_bound import prime_powers


def delay_tails():
    shapes = 0
    for count in (2, 3):
        for z in (1, 2, 4, 8, 16, 32):
            histogram = [0] * (2*z).bit_length()
            for parts in product(range(1, 2*z), repeat=count):
                total = sum(parts)
                if not z <= total < 2*z:
                    continue
                orders = [(a & -a).bit_length()-1 for a in (*parts, total)]
                histogram[max(orders)-min(orders)] += 1
                shapes += 1
            for h in range(len(histogram)):
                # Union over the j slopes and the total, fixing j-1 slopes.
                assert sum(histogram[h:])*2**h <= 2*(count+1)*(2*z)**count
    print("Ordered shapes checked in dyadic delay tails", shapes)


def rows_at_scale(n, scale, powers):
    rows = []
    total = n//scale
    for x in range(1, total//3+1):
        for y in range(x, (total-x)//2+1):
            norm = (x, y, total-x-y)
            a, b, c = (scale*v for v in norm)
            coefficient = comb(n, a)*comb(n-a, b)
            binary = not ((a & b) or (b & c) or (a & c))
            assert binary == bool(coefficient % 2)
            supported = {}
            for p, q in powers:
                zeros = [v//q for v in norm if v % q == 0]
                avoids = bool(coefficient % p)
                if avoids:
                    assert len(zeros) == 2
                    assert sorted(v % q for v in (a, b, c)) == [0, 0, 1]
                if len(zeros) == 2:
                    supported[q] = (zeros, avoids)
            for q, r in combinations(supported, 2):
                part = next(v for v in (a, b, c) if v % (q*r) == 0)
                k = part//(q*r)
                d = (n-1)//(q*r)
                assert k % scale == 0 and 3 <= k < d
                if binary:
                    assert part & ((d-k)*q*r+1) == 0
            if supported:
                rows.append((supported, binary))
    return rows


def split_reconstruction():
    cases = exact = disjoint = exact_binary = 0
    for n in range(3, 1801):
        powers = tuple(prime_powers(n-1))
        if len(powers) < 2:
            continue
        for small in range(3, isqrt(n)+1):
            if n % small or n//small > 90 or gcd(small, n//small) != 1:
                continue
            large = n//small
            left = rows_at_scale(n, small, powers)
            right = rows_at_scale(n, large, powers)
            for u, left_rows, right_rows in ((small, left, right), (large, right, left)):
                for a_powers, binary in left_rows:
                    for c_powers, _ in right_rows:
                        for q, ((x, y), avoids_q) in a_powers.items():
                            for r, (_, avoids_r) in c_powers.items():
                                if q == r:
                                    continue
                                m = u*q
                                d = (n-1)//(q*r)
                                assert u >= 2*r+1 and m > 2*q*r
                                assert 2*(x+y) < d and m*(x+y) <= n-1
                                assert gcd(q, m//q) == 1
                                residue = u*pow(u, -1, q)
                                assert n % m == residue and 2*(n-1) < d*m
                                if binary:
                                    assert m*x & (m*y) == 0
                                cases += 1
                                exact += avoids_q and avoids_r
                                disjoint += binary
                                exact_binary += binary and avoids_q and avoids_r
    assert min(cases, exact, disjoint) > 0
    print("Oriented split controls", cases, "both coefficients avoid their prime", exact,
          "binary first row", disjoint, "both conditions", exact_binary)


def exponent_comparison():
    # Rational upper bound for lambda_3 and lower bound for phi.
    assert 1466**3 > 1466**2*1000+1000**3
    assert 1618**2 < 1618*1000+1000**2
    assert 1466**1000 < 2**552 * 1000**1000
    assert 1618**1000 > 2**694 * 1000**1000
    assert 19*552+3000 < 20*694


if __name__ == "__main__":
    exponent_comparison()
    delay_tails()
    split_reconstruction()
