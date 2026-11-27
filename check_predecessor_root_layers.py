"""Independently test the second-prime residue calculation and the small-root table."""

from itertools import combinations_with_replacement
from math import comb, gcd, isqrt

from check_endpoint_pair_bound import prime_powers
from check_static_power_transfer import prime_powers as powers_through


def eligible(d, m, h, base, top):
    """Test the strongest residue <= 4 exclusion at a full quadratic layer."""
    return (h > 0 and m > 0 and h < top and top > d and d*h <= top
            and d*d*(h*h+4*m) <= 32*top
            and gcd(base, d) == 1 and (d*d-1) % m < 4)


def small_cases(h):
    found = set()
    for m in range(1, h*h//4+1):
        for u in range(1, m+1):
            t, rem = divmod(m, u)
            if rem or t+u != h or t == u or t < 2:
                continue
            for b in range(1, t//2+1):
                if b*u == 1 or gcd(b, t) != 1:
                    continue
                slots = range(b, t, b)
                for i, j in combinations_with_replacement(slots, 2):
                    if i+j > t or (i == j and b % 2 == 0):
                        continue
                    for d in range(m):
                        if h % gcd(d, m):
                            continue
                        k = d*d % m or m
                        if b*k > 2*t or (b*k == 2*t and t*t > h*(i+j)):
                            found.add((t, u, i, j, b, d))
    return found


def check_q_layers():
    # The projection only needs a q-power R and a positive integer P.
    # Enumerate actual upper parts of the second row, including layers below R.
    powers = powers_through(1200)
    seen = passing = lower = equality = 0
    for t in range(2, 9):
        for u in range(1, 9):
            h, m = t+u, t*u
            for R, q in powers.items():
                for d in range(5, 37):
                    P, rem = divmod(d*R-h, m)
                    if rem or P <= 0:
                        continue
                    T, W = 1+t*P, 1+u*P
                    for B in range(1, t//2+1):
                        if T % B:
                            continue
                        for A in range(1, 5):
                            if W % A or gcd(A, B) != 1:
                                continue
                            U, V = B*(W//A), A*(T//B)
                            for z in range(1, (U-1)//(2*R)+1):
                                for w in range(z, (U-1)//R-z+1):
                                    # Choose i=j=B, x=y=A. We use only the
                                    # numerical hypotheses of the projection.
                                    if 2*A*(z+w) >= d:
                                        continue
                                    entries = (V*z, V*w, d*P-V*(z+w))
                                    assert min(entries) >= 0
                                    H = R
                                    while H > d*h and 2*H >= d*(t*t-2*h):
                                        k = (d*d*(R//H)-1) % m+1
                                        assert m*(d*P % H) == k*H-d*h
                                        seen += 1
                                        if sum(a % H for a in entries) == d*P % H:
                                            passing += 1
                                            lower += H < R
                                            equality += B*k == 2*t
                                            assert k >= 5 and B*k >= 2*t
                                            assert B*k > 2*t or t*t > 2*B*h
                                        H //= q
    assert passing and lower
    print(f"Second-prime projection: {seen} layers, {passing} passing additions, "
          f"{lower} passing below the full power, {equality} on the equality boundary.")


def check_grouped():
    cases = 0
    for prime in (3, 5):
        for W in range(prime-1, 2000, prime):
            parts = (prime*W, W+1)
            assert comb(sum(parts), parts[1]) % prime == 0
            cases += 1
    d, P, R = 1001, 1566731, 9391
    assert all(list(prime_powers(v)) == [(v, v)] for v in (P, R))
    m, h = divmod(d*R, P)
    n = 1+d*P*R
    assert n == (1+2*P)*(1+3*P) and eligible(d, m, h, P, P)
    assert isqrt(h*h-4*m)**2 == h*h-4*m
    print(f"Grouped binomials: {cases}. Factored quadratic example: {n}.")


if __name__ == "__main__":
    for h in range(1, 9):
        expected = {(2, 5, 1, 1, 1, 3), (2, 5, 1, 1, 1, 7)} if h == 7 else set()
        assert small_cases(h) == expected
    print("Small remainders 1 through 8: two necessary tuples, both at 7.")
    check_q_layers()
    check_grouped()
