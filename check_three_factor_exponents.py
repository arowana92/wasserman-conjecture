"""Verify the digit counts, shared cut, and explicit progressions with integer arithmetic."""

from itertools import product
from math import comb, gcd


def no_carry(a, b, base):
    while a or b:
        a, da = divmod(a, base)
        b, db = divmod(b, base)
        if da + db >= base:
            return False
    return True


def blocks():
    total = 0
    for p, n in [(2, 7), (3, 4), (5, 3), (7, 2)]:
        for s, t in product(range(1, 12), repeat=2):
            if (s+t) % p == 0:
                continue
            if t % p == 0:
                s, t = t, s
            u, j = s, 0
            while u % p == 0:
                u //= p
                j += 1
            for low in range(p):
                codes = set()
                for block in range(p**n):
                    mult = low + p*block
                    true_digits, code = [], []
                    for i in range(1, n+1):
                        place = p**i
                        a, b = (v*mult//place % p for v in (s, t))
                        if a+b >= p:
                            break
                        digit = u*mult//place % p
                        if j == 0 or i <= j:
                            code.append(digit)
                        else:
                            prefix = mult % place
                            allowed = sorted(u*(prefix+v*place)//place % p
                                             for v in range(p)
                                             if t*(prefix+v*place)//place % p
                                             <= p-1-true_digits[-j])
                            rank = allowed.index(digit)
                            assert rank <= digit and rank+code[-j] < p
                            code.append(rank)
                        true_digits.append(digit)
                    else:
                        word = tuple(code)
                        assert word not in codes
                        codes.add(word)
                if not j:
                    assert len(codes) == ((p+1)//2)**n
                total += 1
    return total


def exact_cuts():
    count = 0
    for d, q, r in product(range(2, 16), (2, 3, 5, 7), (2, 3, 5, 7)):
        if q == r or gcd(d, q*r) != 1:
            continue
        for b, c in product(range(1, 4), repeat=2):
            Q, R = q**b, r**c
            n = 1+d*Q*R
            if n > 5000:
                continue
            for x in range(1, d):
                coefficient = comb(n, x*Q*R)
                assert (coefficient % q != 0) == no_carry(x*R, (d-x)*R, q)
                assert (coefficient % r != 0) == no_carry(x*Q, (d-x)*Q, r)
                count += 1
    return count


def check_progressions():
    count = 0
    for q, r in ((2, 3), (2, 17), (3, 2), (5, 2), (7, 3)):
        stride = 2 if q == 2 else q-1
        for h, first, n in product((1, 2, 6, 19), (1, 7, 31), (1, 2, 3)):
            low = q
            while pow(r, stride*h, low*q) == 1:
                low *= q
            W = low*q**n
            lifts = {pow(r, first+stride*h*i, W) for i in range(q**n)}
            u = (-pow(r, first, low)) % low
            assert 0 < u < low and len(lifts) == q**n and max(lifts) == W-u
            for D in range(2, q**n + 1):
                for x in range(1, D):
                    assert x*(W-u) % W + (D-x)*(W-u) % W >= W
            count += 1
    assert pow(3, 1291, 8192) == 8192-5 and pow(3, 2048, 8192) == 1
    return count


if __name__ == "__main__":
    print("digit-block tests", blocks())
    print("exact binomial cuts", exact_cuts())
    print("progressions", check_progressions())
