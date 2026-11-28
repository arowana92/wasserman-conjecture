"""Check independent Fermat exponents using finite groups and exact cut arithmetic."""

from math import comb, gcd, isqrt


def prime_factors(n):
    factors = []
    for p in range(2, isqrt(n) + 1):
        if n % p == 0:
            factors.append(p)
            while n % p == 0:
                n //= p
    if n > 1:
        factors.append(n)
    return factors


def order(a, m):
    assert gcd(a, m) == 1
    t = a % m
    n = 1
    while t != 1:
        t = t*a % m
        n += 1
    return n


def admissible_at_weight_prime(P, b, exponent, ell, period):
    # Allow independent input units at P and ell. For prime input q all
    # combinations also satisfy the necessary congruences locally.
    target = -pow(pow(2, b, ell)*P, -1, ell) % ell
    images = {pow(v, exponent, ell) for v in range(1, ell)}
    assert pow(2, period, ell) in images
    return target in images


def inspect(m, exponent):
    P = 2**m + 1
    assert prime_factors(P) == [P]
    mask = 4*(exponent & -exponent) - 1
    images = {pow(t, exponent, P) for t in range(1, P)}
    local = []
    for b in range(1, 2*m + 1):
        for t in images:
            x = pow(t*pow(2, b, P), -1, P)
            if x & (P-x) & mask:
                continue
            for ell in prime_factors(x):
                if ell != 2 and admissible_at_weight_prime(P, b, exponent, ell, 2*m):
                    local.append((b, t, x, ell))
    return local


def direct_cuts():
    count = 0
    for P in (5, 17, 257):
        for t in range(1, 12, 2):
            for b in (1, 2, 3):
                N = 1 + P*(2**b)*t
                for x in range(1, P):
                    parity = ((x*t) & ((P-x)*t)) == 0
                    k = x*(2**b)*t
                    assert parity == ((k & (N-k)) == 0)
                    if N <= 1500:
                        assert parity == (comb(N, k) % 2 == 1)
                    count += 1
    return count


def main():
    for m in (8, 16):
        P = 2**m + 1
        E = 2**(m-3)
        x = 2**(m-1) + 1
        assert order(2, P) == 2*m and inspect(m, E) == []
        for ell in prime_factors(x):
            assert order(4, ell) % 2 and (ell-1)//gcd(E, ell-1) % 2
        print(f"P={P}, E={E}: no admissible cut")

    rows = inspect(16, 4096)
    assert {b for b, _, _, _ in rows} == set(range(4, 33, 4))
    assert {(x, ell) for _, _, x, ell in rows} == {(49153, 13)}
    assert all(pow(t, 8, 65537) == 65536 for _, t, _, _ in rows)
    assert 49153 == 13*19*199 and 32769 == 3**2*11*331
    # Enumerate true low-bit powers independently of the group calculations.
    low = {(q % 8, pow(q, 4096, 32768)) for q in range(1, 32768, 2)}
    assert {r for r, t in low if (49153*t & 16384*t & 32767) == 0} == {3, 5}
    print(f"E=4096: {len(rows)} local pairs (b mod 32, power residue) at x=49153, ell=13")
    print(f"{direct_cuts()} binary cut comparisons")


if __name__ == "__main__":
    main()
