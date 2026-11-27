"""Verify paired Legendre layers with signed modular representatives."""

from itertools import product
from math import comb, factorial


def mass(d, u):
    if u % d:
        return 1
    return 2 if u < 0 else 0


def layer(row, m):
    return sum(row)//m - sum(x//m for x in row)


def check_local():
    counts = {2: 0, 3: 0}
    for p in (2, 3, 5, 7):
        for s in (1, 2, 3):
            B = p**s
            for d in range(2, B+2):
                if (B+1) % d:
                    continue
                for j in (1, 2, 3):
                    M = p**j
                    modulus = M*B
                    inverse = pow(d, -1, modulus)
                    values = sorted({v for v in (-M+1, -d, -2, -1, 0, 1, 2, d, M-1)
                                     if abs(v) < M})
                    for m in (2, 3):
                        for us in product(values, repeat=m):
                            if abs(sum(us)) >= M:
                                continue
                            row = [(u*inverse) % modulus + i*modulus for i, u in enumerate(us)]
                            measured = layer(row, M) + layer(row, modulus)
                            predicted = sum(mass(d,u) for u in us) - mass(d,sum(us))
                            assert measured == predicted
                            counts[m] += 1
    print("Paired-layer identities:", counts)


def check_factorials():
    count = 0
    for n in range(3, 91):
        for a in range(1, n//3+1):
            for b in range(a, (n-a)//2+1):
                row = (a, b, n-a-b)
                coefficient = factorial(n)
                for x in row:
                    coefficient //= factorial(x)
                for p in (2, 3, 5):
                    power, valuation = p, 0
                    while power <= n:
                        valuation += layer(row, power)
                        power *= p
                    assert coefficient % p**valuation == 0
                    assert coefficient % p**(valuation+1)
                count += 1
    print("Exact trinomial comparisons:", count)


def check_exponents():
    for ell, s in ((11, 5), (19, 9), (43, 7)):
        assert pow(2, s, ell) == ell-1
        for f in ((1 << (s+1))-1, 3*(1 << (s+1))-1):
            Q = ell**f
            assert Q % 8 == 3 and pow(ell, f+1, 1 << (s+3)) == 1
            for a in range(1, 5):
                cut = (a*Q, (5-a)*Q)
                assert layer(cut, 8) + layer(cut, 1 << (s+3)) == 1
            assert layer((2*Q, 3*Q), 1 << (s+3)) == 1
    # A small full block has exactly half its layers active.
    p, s, d, j = 2, 3, 9, 3
    for v in (-1, 1):
        modulus = p**(j+2*s+1)
        Q = v * pow(d, -1, modulus) % modulus
        cut = (2*Q, 3*Q)
        assert sum(layer(cut,p**i) for i in range(j,j+2*s)) == s
        assert comb(5*Q,2*Q) % p**s == 0
    print("Exponent progressions and full blocks verified")


if __name__ == "__main__":
    check_local()
    check_factorials()
    check_exponents()
