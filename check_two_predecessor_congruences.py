"""Check the reconstruction and constants in two_predecessor_congruences.md."""
from fractions import Fraction as F
from itertools import combinations, product
from math import comb, gcd, prod

from check_binary_windows import LN2_HI, ln_bounds, log2_upper, tangent
from check_endpoint_pair_bound import prime_powers, rows


def check_exponents():
    for a, numbers in (
        (3, (86884, 72186, 64567, 72802, 86890)),
        (2, (93743, 62837, 61750, 64982, 93750)),
    ):
        k, lam, r0, s0, target = (F(n, 100000) for n in numbers)
        assert max(k, 2*(1-k)) < target
        for r, s in product((0, 1), repeat=2):
            edge = (2-lam*(1+F(1, a)))*tangent(r, r0) + 2*lam*tangent(F(r, a), r0/a)
            mid = (1-lam)*tangent(s, s0) + lam*tangent(F(s, a), s0/a)
            assert (1-k)*edge + (2*k-1)*mid < target

    theta = F(29, 40)
    k = theta-F(3, 10000)
    cost = ((1-theta)*tangent(F(1, 3), F(1, 3))
            + theta*(tangent(F(1, 6), F(1, 6))+log2_upper(F(5))/7500)
            + F(2, 3)*(theta-k))
    tail = theta*(1-2/(2500**2*LN2_HI))
    assert max(k, cost, tail) < theta-F(1, 10**7)

    rho, saving = F(1, 2000), F(1, 10**8)
    cases = (
        (1, 2, 2, 0, F(3, 8), F(1, 1000)),
        (1, 3, 2, 0, F(1, 2), F(1, 1000)),
        (1, 2, 2, F(4, 25), F(4, 25), F(1, 1000)),
        (2, 3, 2, F(1, 10000), F(1, 256), F(1, 10000)),
        (2, 3, 3, F(1, 10000), F(1, 256), F(1, 10000)),
    )
    for j, t, p, beta, gamma, eps in cases:
        b = (1-rho-gamma)/j
        mu = (p-1)*(F(1, 2)+beta+eps)/t
        if p == 2:
            cost = tangent(mu, mu)
        else:
            u = F(29, 100)
            cost = (ln_bounds(1+u+u*u)[1] + mu*ln_bounds(1/u)[1])/ln_bounds(F(3))[0]
        assert cost + (1-F(1, t))*(1-b) < 1-saving
        assert 2*eps**2/ln_bounds(F(p))[1] > saving
    print("Global, prescribed-bit, pair and three-row exponent bounds verified")


def recover(bins, v, n):
    """Use two congruences, a factor of the recovered product, then the third."""
    (d1, z1), (d2, z2), (d3, z3) = bins
    w = v*d3
    r1 = pow(z1*d2, -1, d1)
    r2 = pow(z2*d1, -1, d2)
    r = r1 + d1*((r2-r1)*pow(d1, -1, d2) % d2)
    assert w % (d1*d2) == r
    assert w//v == d3
    r3 = pow(v*d1*d2, -1, d3)
    assert z3 % d3 == r3
    total = d1*d2*d3
    restored = w*(z1*d2+z2*d1) + v*z3*d1*d2
    assert restored == n and (n-1) % total == 0
    # U/T controls the number of possibilities in the third residue class.
    if n//v <= 2*total:
        assert z3 in (r3, r3+d3, r3+2*d3)
    return d1 > 1


def coefficients(limit=120):
    count = partial = nonempty = 0
    for n in range(3, limit+1):
        weight = tuple(prime_powers(n))
        pred = tuple(prime_powers(n-1))
        for row, value in rows(n):
            v = prod(w//gcd(w, value) for p, w in weight)
            bins = [1, 1, 1]
            for p, w in pred:
                retained = w//gcd(w, value)
                if retained > 1:
                    residues = [a % retained for a in row]
                    assert sorted(residues) == [0, 0, 1]
                    bins[residues.index(1)] *= retained
                    partial += retained < w
            t = prod(bins)
            assert n*(n-1)//(v*t) == gcd(n*(n-1), value)
            assert all(a % (v*(t//d)) == 0 for a, d in zip(row, bins))
            parameters = sorted((d, a//(v*(t//d))) for a, d in zip(row, bins))
            nonempty += recover(parameters, v, n)
            count += 1
    print(f"{count:,} rows tested with exact coefficients, {partial:,} partially "
          f"retained predecessor powers and {nonempty:,} rows with three nonempty bins")


def large_exact_row():
    row = (86617810, 45491401, 90679143)
    n = sum(row)
    v = t = 1
    bins = [1, 1, 1]
    for m in (n, n-1):
        for p, w in prime_powers(m):
            level, valuation = p, 0
            while level <= n:
                valuation += n//level - sum(z//level for z in row)
                level *= p
            retained = w//gcd(w, p**valuation)
            if m == n:
                v *= retained
            elif retained > 1:
                bins[next(i for i, z in enumerate(row) if z % retained == 1)] *= retained
                t *= retained
    assert n == 222788354 and n*(n-1)//(v*t) < n
    parameters = sorted((d, z//(v*(t//d))) for z, d in zip(row, bins))
    assert recover(parameters, v, n)
    print(f"Reconstructed the three nonempty bins of the exact row at {n:,}")


def constructed():
    tested = small = 0
    for ds in combinations(range(2, 25), 3):
        t = prod(ds)
        if any(gcd(a, b) != 1 for a, b in combinations(ds, 2)):
            continue
        for v in range(1, 25):
            if gcd(v, t) != 1:
                continue
            zs = [pow(v*(t//d), -1, d) for d in ds]
            u = sum(z*(t//d) for d, z in zip(ds, zs))
            assert recover(tuple(zip(ds, zs)), v, u*v)
            tested += 1
            small += u <= t
    examples = 0
    for p in (5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41):
        for m in range(2, p):
            n = m*p+1
            value = n*comb(n-1, p)
            assert gcd(n*(n-1), value) <= n*m and value % p != 0
            examples += 1
    print(f"{tested:,} constructed three-bin reconstructions, {small:,} with U <= T, "
          f"and {examples:,} exact upper-bound examples")


if __name__ == "__main__":
    check_exponents()
    coefficients()
    large_exact_row()
    constructed()
