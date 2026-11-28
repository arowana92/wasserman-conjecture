"""Exercise each low-residue implication using independent base-p arithmetic."""

from itertools import combinations_with_replacement, product
from math import gcd


def fits(part, total, p):
    while part:
        if part % p > total % p:
            return False
        part //= p
        total //= p
    return True


def clear(parts, p):
    power = p
    while power <= p * sum(parts):
        if sum(t % power for t in parts) >= power:
            return False
        power *= p
    return True


def main():
    pred = weight = companion = 0
    for p in (2, 3, 5, 7):
        for m in range(1, 7):
            a = p**m
            for P in (2, 3, 4, 5, 7, 8, 9, 11, 13, 17, 19, 23, 25, 27, 31):
                if P % p == 0:
                    continue
                for t in range(1, 161):
                    if t % p == 0:
                        continue
                    d = P*t % a
                    if P*d > a:
                        continue
                    for x in range(1, P):
                        if gcd(x, P) != 1 or not fits(x*t % a, d, p):
                            continue
                        assert x*d == P*(x*t % a) and d % P == 0
                        pred += 1

                for R in (p, p*p):
                    modulus = R*a
                    for d in range(1, min(a//P, 100) + 1):
                        if d % p == 0:
                            continue
                        # All W in this class give N = R*d modulo R*a.
                        W = (R*d-1)*pow(P, -1, modulus) % modulus
                        for x in range(1, P):
                            if gcd(x, P) != 1 or not fits(x*W % modulus, R*d, p):
                                continue
                            e = x*W % modulus // R
                            assert 0 < e < d and (R*d-1)*x == R*P*e
                            assert (R*d-1) % P == 0
                            weight += 1

            for u in range(3, 29):
                if gcd(u, p) != 1:
                    continue
                for d in range(1, min(a//(u+1), 60) + 1):
                    if d % p == 0:
                        continue
                    z = d*pow(u, -1, a) % a
                    for g, h in combinations_with_replacement(range(1, u), 2):
                        k = u - g - h
                        if k < h:
                            continue
                        parts = (g, h, k)
                        e = tuple(z*x % a for x in parts)
                        if sum(e) != d or not clear(e, p):
                            continue
                        assert all(x*d == u*y and y > 0 for x, y in zip(parts, e))
                        assert gcd(d, u) >= 3
                        companion += 1
    assert min(pred, weight, companion) > 0
    for i, j in product(range(1, 11), (1, 101, 201)):
        N = 1 + 11*5**i*23**j
        assert (N-1)//5**i % 125 == 3
    print(f"surviving local configurations: predecessor {pred}, weight cut {weight}, companion {companion}")


if __name__ == "__main__":
    main()
