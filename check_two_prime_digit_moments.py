"""Compare exact two-prime exponent-period counts with independent floor tests."""
from fractions import Fraction
from itertools import product


def valuation(n, p):
    v = 0
    while n % p == 0:
        v += 1
        n //= p
    return v


def moment(p):
    direct = Fraction(sum((d-m+1)**2 for d in range(p) for m in range(d+1)), p**4)
    closed = Fraction((p+1)**2*(p+2), 12*p**3)
    assert direct == closed
    return closed


def room(u, m, p, length):
    choices = 1
    for _ in range(length):
        choices *= max(0, 1 + u % p - m % p)
        u //= p
        m //= p
    return choices


def concentration(count, modulus):
    return Fraction(((count+modulus-1)//modulus)*modulus, count)


def squared_bound(p, q, a, b, v, w, j, k):
    return (concentration(q**(b-w), p**j) * moment(p)**(j-v)
            * concentration(p**(a-v), q**k) * moment(q)**(k-w))


def exact_counts(p, q, r, step, start, a, b, j, k):
    P, Q = p**a, q**b
    v, w = valuation(r**step-1, p), valuation(r**step-1, q)
    assert 1 <= v <= j <= a and 1 <= w <= k <= b
    assert (p != 2 or v >= 2) and (q != 2 or w >= 2)
    T = p**(a+j-v) * q**(b+k-w)
    x0 = pow(Q * r**start, -1, p**v)
    y0 = pow(P * r**start, -1, q**w)
    total = 0
    for x in range(x0, P, p**v):
        for y in range(y0, Q, q**w):
            rp, rq = pow(x*Q, -1, P), pow(y*P, -1, Q)
            total += (room(Q*rp, y*rp, p, j)
                      * room(P*rq, x*rq, q, k))

    modulus = P * p**j * Q * q**k
    z = pow(r, start, modulus)
    step_factor = pow(r, step, modulus)
    direct = 0
    for _ in range(T):
        x, y = pow(Q*z, -1, P), pow(P*z, -1, Q)
        weight = 1 + P*Q*z
        row = (x*Q*z, y*P*z, weight-x*Q*z-y*P*z)
        # Signed third parts are allowed in this truncated residue test.
        for base, length in ((p, a+j), (q, b+k)):
            if any(sum(n // base**i for n in row) != weight // base**i
                   for i in range(1, length+1)):
                break
        else:
            direct += 1
        z = z * step_factor % modulus
    assert z == pow(r, start, modulus)
    assert total == direct
    assert Fraction(total, T)**2 <= squared_bound(p, q, a, b, v, w, j, k)
    return total, T


def log_bounds(x):
    # Positive-term expansion of log(x) with a geometric remainder bound.
    z = (x-1)/(x+1)
    assert 0 < z < 1
    low = sum(2*z**i/i for i in range(1, 40, 2))
    return low, low + 2*z**41/(41*(1-z*z))


def check_interval_exponents():
    f = Fraction
    log2, log3 = log_bounds(f(2)), log_bounds(f(3))
    first, second = log_bounds(f(8, 3)), log_bounds(f(81, 20))
    k_low = (first[0]/log2[1] + second[0]/log3[1])/2
    k_high = (first[1]/log2[0] + second[1]/log3[0])/2

    # Lambda_2 and Lambda_3 are the roots in (1,2) and (2,3) respectively.
    # For Lambda_3 use 1/(2 sin(pi/14)) = 1 + 2 cos(2 pi/7).
    cubic = lambda x: x**3-2*x*x-x+1
    lo, hi = f(2246979, 10**6), f(2246980, 10**6)
    assert cubic(lo) < 0 < cubic(hi)
    beta_low = log_bounds(lo)[0]/log3[1]
    beta_high = log_bounds(hi)[1]/log3[0]
    phi_upper = f(1618034, 10**6)
    assert phi_upper**2-phi_upper-1 > 0
    assert log_bounds(phi_upper)[1]/log2[0] < beta_low

    u_low = (1-beta_high)/(1+k_high)
    u_high = (1-beta_low)/(1+k_low)
    assert u_low > f(1122, 10000)
    assert 1-2*u_high > f(7755, 10000)
    assert beta_high+u_high < f(8492, 10000)


def main():
    check_interval_exponents()
    checked = 0
    for p, q, r, step, da, db in [
            (2, 3, 5, 2, 2, 2), (2, 5, 3, 4, 1, 1),
            (3, 5, 2, 4, 1, 1), (5, 7, 13, 12, 1, 0)]:
        v, w = valuation(r**step-1, p), valuation(r**step-1, q)
        for a, b in product(range(v, v+da+1), range(w, w+db+1)):
            for j, k in product(range(v, a+1), range(w, b+1)):
                for start in (1, 2):
                    if p**(a+j-v) * q**(b+k-w) > 20000:
                        continue
                    exact_counts(p, q, r, step, start, a, b, j, k)
                    checked += 1
    print(checked, 'joint period counts verified')
    for start in (1, 2):
        print('a=5 b=3, parity', start % 2,
              exact_counts(2, 3, 5, 2, start, 5, 3, 5, 3))
    for a, b, cutoff in [(64, 40, Fraction(713, 10**27)),
                         (128, 81, Fraction(588, 10**53))]:
        assert squared_bound(2, 3, a, b, 3, 1, a, b) < cutoff**2
        print('a, b', a, b, 'failure fraction less than', cutoff)


if __name__ == '__main__':
    main()
