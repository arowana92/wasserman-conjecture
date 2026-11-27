"""Certify entropy exponents and test the digit-sum complement on exact pairs."""

from fractions import Fraction as F
from itertools import combinations_with_replacement
from math import comb, gcd, prod

from check_binary_windows import ln_bounds
from check_endpoint_pair_bound import prime_powers


def digit_bound(p, mean, q):
    entropy = ln_bounds(sum(q**d for d in range(p)))[1] + mean*ln_bounds(1/q)[1]
    return entropy / ln_bounds(F(p))[0]


def constants():
    for p, alpha, q, target in (
        (2, F(0), F(1,3), F(8113,10000)), (3, F(0), F(4343,10000), F(8204,10000)),
        (2, F(1,10), F(3793,10000), F(8561,10000)),
        (3, F(1,10), F(48,100), F(8561,10000)),
        (5, F(1,20), F(5883,10000), F(8561,10000)),
        (7, F(1,40), F(6601,10000), F(8561,10000)),
    ):
        assert digit_bound(p, (p-1)*(1+alpha)/4, q) < target

    # Cut off the density of digits in the weight. Columns are in units of 1/10000.
    for b, p, r0, high, low, bound in (
        (F(15,16), 2, 7296, 3706, 3214, 8511),
        (F(15,16), 3, 7287, 4731, 4215, 8511),
        (F(7,8), 2, 7078, 4128, 3088, 8805),
        (F(7,8), 3, 7058, 5161, 4079, 8805),
        (F(31,48), 2, 6101, 6391, 2553, 9717),
        (F(31,48), 3, 6016, 7337, 3474, 9717),
    ):
        r, u, v, target = (F(n,10000) for n in (r0, high, low, bound))
        assert digit_bound(p, (p-1)*(1-r), u) < target
        assert digit_bound(p, (p-1)*r/3, v) + F(2,3)*(1-b) < target


def digit_sum(n, p):
    count = 0
    while n:
        n, d = divmod(n, p)
        count += d
    return count


def counts():
    pair_count = row_count = 0
    for n in range(3, 37):
        layers = [list(prime_powers(m)) for m in (n, n-1)]
        rows = []
        for a in range(1, n//3+1):
            for b in range(a, (n-a)//2+1):
                parts = (a, b, n-a-b)
                coef = comb(n,a)*comb(n-a,b)
                for p in (2,3,5,7):
                    value = coef
                    v = 0
                    while value%p == 0:
                        value //= p
                        v += 1
                    power = p
                    length = 1
                    while power <= n:
                        power *= p
                        length += 1
                    scores = [digit_sum(z,p) for z in (*parts, power-1-n)]
                    assert sum(scores) == (p-1)*(length+v)
                rows.append((parts, coef))
                row_count += 1
        for ((parts1,c1),(parts2,c2)) in combinations_with_replacement(rows, 2):
            assigned = [[1,1],[1,1]]
            for side, factors in enumerate(layers):
                for _, w in factors:
                    g1, g2 = gcd(c1,w), gcd(c2,w)
                    who = int(g2 < g1)
                    assigned[who][side] *= w//min(g1,g2)
            g = gcd(n*(n-1),gcd(c1,c2))
            assert prod(assigned[0])*prod(assigned[1])*g == n*(n-1)
            for (u,s), parts in zip(assigned, (parts1,parts2)):
                assert all(z%u == 0 and z*(z-1)%s == 0 for z in parts)
            if 1 in parts1:
                assert all(z % (n//gcd(n,c2)) == 0 for z in parts2)
            pair_count += 1
    print(f'Checked {row_count} exact rows and {pair_count} unordered row pairs.')


if __name__ == '__main__':
    constants()
    counts()
