"""Check the quadratic cut, its progression, and the finite-window binary estimate."""

from itertools import product
from math import gcd, lcm

from check_endpoint_pair_bound import prime_powers


def v2(n):
    assert n
    return (n & -n).bit_length()-1


def at_most_phi_power(count, exponent):
    # 2*phi**e = Lucas(e) + sqrt(5)*Fibonacci(e).
    f, g = 0, 1
    for _ in range(exponent):
        f, g = g, f+g
    difference = 2*count-(2*g-f)
    return difference <= 0 or difference*difference <= 5*f*f


def finite_windows():
    blocks = words = positive = 0
    for a, b in product(range(4), repeat=2):
        for c, e in ((1,1), (2,8), (3,15), (12,4), (8,16)):
            for m in range(1, 6):
                size = 1 << m
                mask = size-1
                for prefix in range(size):
                    slope = 2*a*prefix+c
                    other = 2*b*prefix+e
                    a0 = (a*prefix*prefix+c*prefix+7)//size
                    b0 = (b*prefix*prefix+e*prefix+12)//size
                    kept = 0
                    for y in range(size):
                        x = prefix+size*y
                        ah = ((a*x*x+c*x+7)//size) & mask
                        bh = ((b*x*x+e*x+12)//size) & mask
                        assert ah == (a0+slope*y) & mask
                        assert bh == (b0+other*y) & mask
                        kept += (ah & bh) == 0
                    w = max(v2(n) for n in (slope, other, slope+other))
                    assert at_most_phi_power(kept, m+min(w,m))
                    blocks += 1
                    words += size
                    positive += kept
    print(f"Binary tangent windows: {blocks} blocks, {words} upper words, "
          f"{positive} carry-free in the tested window.")


def assigned_groups(row, factorization):
    bins = [1, 1, 1]
    for power in factorization:
        residues = [a % power for a in row]
        if sorted(residues) != [0, 0, 1]:
            return None
        bins[residues.index(1)] *= power
    return bins


def actual_rows():
    seen = legal = unit_bins = even_bins = extra = 0
    # Allow the positive input P to be composite, as the companion lemma
    # needs no primality assumption once the first-row shape is specified.
    for t, u in ((2,5), (3,2), (4,4), (9,6), (11,3), (13,5)):
        M, h = t*u, t+u
        for P in range(2, 32):
            T, W, H = 1+t*P, 1+u*P, M*P+h
            N = T*W
            for B in range(1, t//2+1):
                if T % B:
                    continue
                for S in range(1, t*t//4+1):
                    if H % S or gcd(S, H//S) != 1:
                        continue
                    # Include S if it can divide the two bounded multipliers' lcm.
                    if not any(lcm(i,j) % S == 0 for i in range(B,t,B)
                               for j in range(i,t-i+1,B)):
                        continue
                    powers = [w for _, w in prime_powers(H//S)]
                    for a in range(1, B*W//3+1):
                        for b in range(a, (B*W-a)//2+1):
                            row = [T//B*x for x in (a,b,B*W-a-b)]
                            seen += 1
                            bins = assigned_groups(row, powers)
                            if bins is None:
                                continue
                            d = min(bins)
                            Y = row[bins.index(d)]
                            n, rem = divmod(Y, (T//B)*(H//(S*d)))
                            assert not rem and gcd(n,d) == 1 and H//S >= d**3
                            L, Q = B*S*d, B*S*d*d
                            assert 0 < t*n < L
                            assert (t*n*(H//(S*d))+B*u) % d == 0
                            assert (Y*L, N-Y) == (T*H*n, N-T*H*n//L)
                            E = L-t*n
                            changes = (d*n*(2*t*H-t*t), d*(E*(2*H-h)-M*n))
                            leading = (n*t*M*Q*d, E*M*Q*d)
                            for step in (1, 3):
                                pp = P+step*Q
                                hh, tt = M*pp+h, 1+t*pp
                                yy = tt*hh*n//L
                                nn = 1+pp*hh
                                assert (yy-Y, (nn-yy)-(N-Y)) == tuple(
                                    step*slope+step**2*quad
                                    for slope,quad in zip(changes,leading))
                            legal += 1
                            unit_bins += d == 1
                            even_bins += d % 2 == 0
                            extra += B > 1 or S > 1
    assert min(legal, unit_bins, even_bins, extra) > 0
    print(f"Smallest group: {legal} valid allocations among {seen} actual scaled rows, "
          f"{unit_bins} with group one, {even_bins} with even smallest group, "
          f"{extra} with B or S greater than one.")


if __name__ == "__main__":
    finite_windows()
    actual_rows()
