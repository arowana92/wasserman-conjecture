"""Verify the inverse-residue theorems at all exponents (requires python-flint)."""
from collections import Counter
from functools import lru_cache
from math import gcd, lcm, prod

from flint import fmpz
from check_three_factors import valuation
from prime_certificates import primes_through


@lru_cache(maxsize=32768)
def factor(n):
    result = tuple((int(q), int(e)) for q, e in fmpz(n).factor())
    assert prod(q**e for q, e in result) == n
    assert all(e > 0 and fmpz(q).is_prime() for q, e in result)
    assert len({q for q, _ in result}) == len(result)
    return result


def order(p, d):
    phi = prod((q-1)*q**(e-1) for q, e in factor(d))
    m = phi
    for q, _ in factor(phi):
        while m % q == 0 and pow(p, m//q, d) == 1:
            m //= q
    assert pow(p, m, d) == 1 and all(pow(p, m//q, d) != 1 for q, _ in factor(m))
    return m


def digits_fit(r, t, p):
    if p == 2:
        return r & t == r and bin(r).count("1") >= 4
    total = 0
    while r:
        r, a = divmod(r, p)
        t, b = divmod(t, p)
        if a > b:
            return False
        total += a
    return total >= 4


def carry_free(x, y, p):
    if p == 2:
        return x & y == 0
    while x and y:
        x, a = divmod(x, p)
        y, b = divmod(y, p)
        if a+b >= p:
            return False
    return True


def empty_level(q, p, t, parts, g, scale_base, start, step):
    """All powers p^j, j=start modulo step, are among the surviving lifts."""
    x, y, z = parts
    survivors, w = [0], 1
    for _ in range(20):
        old, w = w, w*q
        if pow(p, step, w) == 1:
            candidates = (pow(p, start, w),)
        elif len(survivors)*q <= 30000:
            candidates = (a+old*i for a in survivors for i in range(q))
        else:
            return None
        survivors = [a for a in candidates if
                     (x*a) % w+(y*a) % w+(z*a+1) % w == (t*a+1) % w
                     and not (q == scale_base and (t*a+1) % w == 0 and g % w)]
        if not survivors:
            return w
    return None


def check_scales(p, t, m, h, factors, parts, counts):
    x, y, z = parts
    g = gcd(x, y)
    choices = [q for q, _ in factor(g) if q != 2]
    if g % 4 == 0:
        choices.append(4)
    for v in choices:
        if gcd(v, p*t) > 1:
            continue
        step = lcm(m, order(p, v))
        for start in range((-h) % m, step, m):
            if (t*pow(p, start, v)+1) % v:
                continue
            counts["progressions"] += 1
            base = 2 if v == 4 else v
            for q in (base, factors[0][0], factors[1][0]):
                if empty_level(q, p, t, parts, g, base, start, step) is not None:
                    break
            else:
                raise AssertionError(("unexcluded progression", p, t, m, h, v, start, step))


def check_cofactor(p, d, counts):
    m = order(p, d)
    t = (p**m-1)//d
    counts["base and quotient pairs"] += 1
    if t <= 1:
        return
    power = p**(m-1)
    for h in range(m-1, -1, -1):
        if power < t:
            break
        r = power % t
        power //= p
        if not digits_fit(r, t, p):
            continue
        counts["digit candidates"] += 1
        factors = factor(t)
        if len(factors) != 2:
            continue
        P, Q = (q**e for q, e in factors)
        x = Q*(r*pow(Q, -1, P) % P)
        y = P*(r*pow(P, -1, Q) % Q)
        if x+y != r or not carry_free(x, y, p) or gcd(x, y) < 3:
            continue
        counts["splits"] += 1
        check_scales(p, t, m, h, factors, (x, y, t-r), counts)


def small_crt():
    checked = passing = 0
    for p in (2, 3, 5, 7):
        for t in range(6, 200):
            fs = factor(t)
            if len(fs) != 2 or t % p == 0:
                continue
            P, Q = (q**e for q, e in fs)
            for j in range(1, 7):
                w = p**j
                n = w*t+1
                # Build the three actual CRT idempotents, without using r.
                row = tuple(((n-1)//v)*pow((n-1)//v, -1, v) for v in (P, Q, w))
                x, y = row[0]//w, row[1]//w
                r = pow(w, -1, t)
                assert (sum(row) == n) == (x+y == r)
                assert x == Q*(r*pow(Q, -1, P) % P)
                assert y == P*(r*pow(P, -1, Q) % Q)
                if sum(row) == n:
                    assert row[2] == w*(t-r)+1
                    predicted = carry_free(x, y, p) and carry_free(r, t-r, p)
                    assert predicted == (valuation(n, row, p) == 0)
                    if predicted:
                        assert digits_fit(r, t, p)
                        passing += 1
                checked += 1
    assert passing > 0
    print(f"{checked:,} direct CRT row checks, including {passing:,} carry-free controls")


if __name__ == "__main__":
    small_crt()
    assert fmpz(2**89-1).is_prime() and fmpz(179).is_prime() and (2**89+1) % 179 == 0
    counts = Counter()
    limits = {2: 30000, 3: 12000, 5: 12000, 7: 12000}
    for p in primes_through(2500):
        for d in range(p+1, limits.get(p, 2500)+1):
            if d % p:
                check_cofactor(p, d, counts)
    print(dict(counts))
