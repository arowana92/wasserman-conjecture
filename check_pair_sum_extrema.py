"""Verify the two-part sum bound with exact coefficients and retained prime powers."""

from collections import Counter
from itertools import combinations
from math import gcd, prod

from check_endpoint_pair_bound import allocate, prime_powers, rows


def bins(parts, modulus):
    result = tuple(gcd(modulus, parts[i]-1, parts[(i+1)%3], parts[(i+2)%3])
                   for i in range(3))
    assert prod(result) == modulus
    return result


def exact_comparisons(bound):
    counts = Counter()
    for n in range(3, bound+1):
        table = list(rows(n))
        weight, predecessor = tuple(prime_powers(n)), tuple(prime_powers(n-1))
        for row, ka in table:
            for companion, kc in table:
                D, J = gcd(n, gcd(ka, kc)), gcd(n-1, gcd(ka, kc))
                u, v, partial_n = allocate(weight, ka, kc)
                s, t, partial_m = allocate(predecessor, ka, kc)
                first, second = bins(row, s), bins(companion, t)
                co = [y//(v*(t//m)) for y, m in zip(companion, second)]
                ell = prod(co)
                assert all(y == v*f*(t//m) for y, f, m in zip(companion, co, second))
                assert (D*u)**3 >= 27*ell*t*t
                for i, k in combinations(range(3), 2):
                    x, y = row[i], row[k]
                    h = x+y
                    d, e, w = first[i], first[k], first[3-i-k]
                    a, b = x//(u*w*e), y//(u*w*d)
                    c = (h-1)//(d*e)
                    assert (x, y, h-1) == (u*w*a*e, u*w*b*d, c*d*e)
                    penalty = w*w*(a*b)**3*c*ell**2
                    lhs = D**6*J**4*h**6*(h-1)
                    rhs = 6**6*(n-1)**4*penalty
                    assert lhs >= rhs
                    # The ratio factors exactly into a two-term and a three-term AM–GM ratio.
                    assert rhs*h**6*(D*u)**6 == lhs*(4*x*y)**3*(27*ell*t*t)**2
                    counts['selected pairs'] += 1
                    counts['partly retained'] += partial_n or partial_m
                    m = n-1
                    # Remove full powers until at most four distinct primes remain.
                    atoms = sorted(power for _, power in predecessor)
                    r = prod(atoms[:max(0, len(atoms)-4)])
                    assert D**3*J**2*r*r*h*h*(h-1) > 4*m*m
                    if (D*J)%2:
                        if n%2 == 0:
                            assert penalty >= (512 if h%2 else 32)
                            if v%2 == 0:
                                assert ell%8 == 0
                                assert (D*u)**3 > 32*ell*t*t
                            else:
                                assert (a*b)%2 == 0 and ell%2 == 0
                            counts['even weights, odd gcd'] += 1
                        elif h%2 == 0:
                            assert penalty >= 4
                            counts['odd weights, even sum, odd gcd'] += 1
    print(dict(counts))


def two_leading_units(bound):
    comparisons = 0
    for n in range(3, bound+1):
        for p, _ in prime_powers(n*(n-1)):
            leading = p
            while leading*p <= n:
                leading *= p
            offset = n-leading
            if not 0 < offset < leading:
                continue
            middle = 1
            while middle*p <= offset:
                middle *= p
            tail = offset-middle
            if tail >= middle:
                continue
            for (x,y,z), value in rows(n):
                if value%p == 0:
                    continue
                assert leading <= z and x*y <= middle*tail
                comparisons += 1
    assert comparisons > 0
    print('Two leading units, exact avoiding rows:', comparisons)
    for i, a in enumerate((5291, 6968, 11193, 15257, 23875, 40331, 51368)):
        target = 6**6*2**(4*i)
        b = 1
        while (2*b)**7 <= target:
            b *= 2
        assert (1000*b)**7 < 1000**7*target < a**7
        assert a < 2000*b and a < 1000*2**(3+i)
        assert (b*(a-1000*b))**3*a <= 729*2**(4*i)*1000**4
        assert (b*(a+1-1000*b))**3*(a+1) > 729*2**(4*i)*1000**4
        if i == 5:
            assert a**7 > 1190**7*target


if __name__ == '__main__':
    exact_comparisons(40)
    two_leading_units(150)
