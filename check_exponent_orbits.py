"""Check the exponent congruences, class counts, and translated-cube arithmetic."""
from collections import Counter
from itertools import product
from math import comb, gcd


def carry_free(x, y, p):
    while x or y:
        if x % p + y % p >= p:
            return False
        x //= p
        y //= p
    return True


def cycle(a, modulus):
    values = [1]
    x = a % modulus
    while x != 1:
        values.append(x)
        x = x*a % modulus
    return values


def check_orbits():
    comparisons = 0
    for p, q, r in ((2, 3, 5), (3, 2, 7), (5, 2, 3), (5, 3, 7)):
        for a in range(1, 4):
            P = p**a
            for h, k in ((1, 1), (6, 2), (P*(p-1), 1), (p-1, P*(p-1))):
                qs, rs = cycle(pow(q, h, P), P), cycle(pow(r, k, P), P)
                fibers = Counter((v*w) % P for v in qs for w in rs)
                assert set(fibers.values()) == {len(qs)*len(rs)//len(fibers)}
                for d in ((P, 2*P) if min(p, q, r) > 2 else (P,)):
                    for B, H in ((2, 3), (11, 7), (5, 21)):
                        reached = set()
                        for (i, u), (j, v) in product(enumerate(qs), enumerate(rs)):
                            x0 = pow(q*r*u*v, -1, P)
                            for x in range(x0, d, P):
                                ib = [m for m in range(i, B, len(qs))
                                      if carry_free(x*q**(1+h*m), (d-x)*q**(1+h*m), r)]
                                jc = [m for m in range(j, H, len(rs))
                                      if carry_free(x*r**(1+k*m), (d-x)*r**(1+k*m), q)]
                                reached.update((m, n, x) for m in ib for n in jc)
                        direct = set()
                        for i, j in product(range(B), range(H)):
                            Q, R = q**(1+h*i), r**(1+k*j)
                            x0 = pow(Q*R, -1, P)
                            for x in range(x0, d, P):
                                s, t = x*Q*R, 1+(d-x)*Q*R
                                if carry_free(s, t, q) and carry_free(s, t, r):
                                    direct.add((i, j, x))
                                if s+t <= 4000:
                                    coefficient = comb(s+t, s)
                                    assert (coefficient % q != 0 and coefficient % r != 0) == (
                                            (i, j, x) in direct)
                        assert reached == direct
                        comparisons += 1
    return comparisons


def check_final_valuations():
    cases = [(1, 4, 1, 1), (1, 24, 1, 6), (1, 4, 2, 1), (1, 6, 2, 2),
             (2, 1, 3, 1), (2, 7, 5, 2), (6, 4, 4, 2)]
    for a, b, c, expected in cases:
        d, w = 5**a, 2**b*3**c
        n, part = 1+d*w, pow(w, -1, d)*w
        for p, target in ((2, 0), (3, 0), (5, expected)):
            place, val = p, 0
            while place <= n:
                val += n//place-part//place-(n-part)//place
                place *= p
            assert val == target
    return len(cases)


def check_farey_layers():
    checks = 0
    for d in range(2, 55):
        for m in (16, 27, 49, d, 6*d):
            for n in range(m):
                low = d*n % m
                direct = all(x*n % m > low for x in range(1, d))
                u = d*n//m
                formula = gcd(u, d) == 1 and low * (-pow(u, -1, d) % d) < m
                assert direct == formula
                checks += 1
    d, M = 5**10, 2**40
    assert 2*sum(M//(d*t) for t in range(1, min(d-1, M//d)+1)
                 if gcd(d, t) == 1) == 2195614
    return checks


def lattice_basis():
    # Compute log_2(3) in the smooth-order unit group modulo 5^1024.
    m = 5**1023
    M = 5*m
    ell, order, modulus = 3, 4, 5
    z, step = 8*pow(3, -1, M) % M, 16
    for _ in range(1, 1024):
        modulus *= 5
        for digit in range(5):
            if z % modulus == 1:
                break
            z = z*step % M
        else:
            raise ArithmeticError('missing discrete logarithm digit')
        ell += digit*order
        order *= 5
        step = pow(step, 5, M)
    assert pow(2, ell, M) == 3
    U, V = 4*3**31, 2**62
    slope = -(V//4)*ell*pow(U//4, -1, m) % m
    u, v = (U*m, 0), (U*slope, V)
    dot = lambda x, y: x[0]*y[0] + x[1]*y[1]
    while True:
        if dot(u, u) > dot(v, v):
            u, v = v, u
        n = (2*dot(u, v)+dot(u, u))//(2*dot(u, u))
        if n == 0:
            break
        v = v[0]-n*u[0], v[1]-n*u[1]
    options = []
    for i, j in product(range(-9, 10), repeat=2):
        x = i*u[0]+j*v[0], i*u[1]+j*v[1]
        if min(x) > 0:
            options.append((x, i, j))
    _, x, y = min((max(*x, *y), x, y) for x, i, j in options for y, k, l in options
                  if abs(i*l-j*k) == 1)
    assert abs(x[0]*y[1]-x[1]*y[0]) == U*V*m
    for h, k in (x, y):
        assert h % (2*3**31) == k % V == 0
        assert pow(2, h, M)*pow(3, k, M) % M == 1
    assert max(x[0]+y[0], x[1]+y[1])+1024 < 10**375
    return x, y


def check_transfers():
    # Check the exact numerator identity using small, nonzero shifts of all axes.
    checks = 0
    for a, b, c in product(range(1, 4), repeat=3):
        W, D = 2**b*3**c, 5**a
        d = W * 2**2 * 3**2
        A = len(cycle(5, d))
        new_D = 5**(a+A)
        x, new_x = pow(W, -1, D), pow(W, -1, new_D)
        t = (W*x-1)//D
        assert W*new_x == 1+t*new_D and (new_x-x) % 36 == 0
        modulus = new_D*25
        # Full orders are unnecessary to exhibit one nonzero period.
        B = 4*(modulus//5)*6
        C = 4*(modulus//5)*2
        new_W = W*pow(2, B, modulus)*pow(3, C, modulus) % modulus
        assert pow(new_W, -1, new_D) == new_x
        assert (new_x*new_W-1)//new_D % 25 == t % 25
        checks += 1
    A0 = 2**1086 * 3**1055
    assert pow(5, A0, 2**1088 * 3**1056) == 1
    return checks


if __name__ == '__main__':
    print(check_orbits(), 'exact rectangle comparisons')
    print(check_transfers(), 'translated triples checked')
    print(check_final_valuations(), 'binomial valuation triples checked')
    print(check_farey_layers(), 'rational carry layers checked')
    for vector in lattice_basis():
        print('positive period vector', *vector)
