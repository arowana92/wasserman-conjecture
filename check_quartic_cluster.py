"""Verify the carry and resultant arguments for clusters of prime powers."""
from collections import Counter
from itertools import combinations, product
from math import gcd, lcm, prod

from flint import fmpz_poly

from check_scaled_row_blocks import certify_prime
from check_translated_cubic_bounds import resultant
from search_endpoint import atoms, smallest_factors


def middle_layer():
    powers = []
    sieve = smallest_factors(50)
    for w in range(2, 51):
        factors = atoms(w, sieve)
        if len(factors) == 1:
            powers.append((factors[0][0], w))
    counts = Counter()
    for (p, x), (q, y), (r, z) in combinations(powers, 3):
        if len({p, q, r}) < 3:
            continue
        block = x*y*z
        for k1, k2, k3 in product(range(1, 5), repeat=3):
            normalized = (k1*y*z, k2*x*z, k3*x*y)
            u = sum(normalized)
            if gcd(u, block) != 1:
                continue
            v = pow(u, -1, block)
            row = tuple(v*a for a in normalized)
            n = sum(row)
            m = (n-1)//block
            budget = (m*x*z) % y
            first = (row[0]//y) % y
            assert (k1 + k2*(y-x)*first) % y == 0
            # Check the actual factorial-floor term at the next block boundary.
            if n//(y*y) == sum(a//(y*y) for a in row):
                counts["no carry at this layer"] += 1
                assert first <= budget and y <= k1+k2*(y-x)*budget
            else:
                counts["carry at this layer"] += 1
            if y > k1+k2*(y-x)*budget:
                assert n//(y*y) > sum(a//(y*y) for a in row)
                counts["strict exclusions"] += 1
    assert min(counts.values()) > 0 and len(counts) == 3
    return dict(counts)


def polynomial_nodes():
    counts = Counter()
    x = fmpz_poly([0, 1])
    # Check all individual coefficient splits, including the square F+1.
    for shifts in [(0, 1, 2, 3), (0, 2, 3, 5), (0, 1, 3, 6)]:
        h = max(shifts)
        f = prod((x+s for s in shifts), start=fmpz_poly([1]))
        for offset in (0, 1):
            total = f + offset
            coefficients = [int(t) for t in total.coeffs()]
            for values in product(*(range(t+1) for t in coefficients)):
                part = fmpz_poly(list(values))
                if part.degree() <= 0:
                    continue
                roots = [s for s in shifts if part(-s) == 0]
                phi = prod((x+s for s in roots), start=fmpz_poly([1]))
                quotient, rem = divmod(part, phi)
                assert not rem
                missing = len(shifts)-len(roots)
                for s in shifts:
                    val = int(part(-s))
                    bound = s**4+coefficients[2]*s*s
                    assert -bound <= val <= bound+offset
                    if val in (0, 1) and any(part(-t) in (0, 1) and abs(s-t) >= 2
                                              and part(-t) != val for t in shifts):
                        raise AssertionError("distinct unit components at separated nodes")
                if not roots:
                    continue
                target = f+1 if offset else f-1
                res = resultant(quotient.coeffs(), target.coeffs())
                if res:
                    assert res <= (2*3**4)**4 * (h+1)**(4*missing)
                    counts["nonzero resultants"] += 1
                else:
                    pairs = target.factor()[1]
                    assert len(pairs) == 1 and pairs[0][1] == 2
                    factor = pairs[0][0]
                    reduced, rem = divmod(quotient, factor)
                    assert not rem
                    res = resultant(reduced.coeffs(), factor.coeffs())
                    assert res and res*res <= (2*3**4)**4 * (h+1)**(4*missing)
                    counts["square-factor divisions"] += 1
    return dict(counts)


def linear_bins():
    counts = Counter()
    x = fmpz_poly([0, 1])
    h = 6
    for s in (2, 4, 6):
        t, u = [v for v in (2, 4, 6) if v != s]
        vpoly = (x+t)*(x+u)
        phi = x*(x+s)
        f = phi*vpoly
        # Actual normalized small part polynomials with a nonnegative product
        # by phi. Allow the two parts to coincide.
        parts = [a*x+b for a in range(6) for b in range(15) if a or b]
        for i, first in enumerate(parts):
            for second in parts[i:]:
                third = vpoly-first-second
                if any(v < 0 for v in (phi*third).coeffs()):
                    continue
                a, b = int(first[1]), int(first[0])
                g, e = int(second[1]), int(second[0])
                p = 10**24+37
                m = int(f(p)-1)
                actual = [int(first(p)), int(second(p)), int(third(p))]
                bins = [gcd(m, gcd(actual[j], actual[k])) for j, k in combinations(range(3), 2)]
                assert lcm(*bins) < (8*h)**11
                if not a or not g:
                    counts["constant"] += 1
                elif a*e == b*g:
                    counts["proportional"] += 1
                else:
                    d1 = b*b-(t+u-g)*a*b+(t*u-e)*a*a
                    d2 = e*e-(t+u-a)*g*e+(t*u-b)*g*g
                    assert abs(d1) == resultant(first.coeffs(), (vpoly-second).coeffs())
                    assert abs(d2) == resultant(second.coeffs(), (vpoly-first).coeffs())
                    assert max(abs(d1), abs(d2)) <= 72*h**4
                    for poly, quad, result in ((first, vpoly-second, d1), (second, vpoly-first, d2)):
                        if result:
                            assert result % gcd(int(poly(p)), int(quad(p))) == 0
                        else:
                            root, rem = divmod(-int(poly[0]), int(poly[1]))
                            assert not rem and abs(root) <= 4*h
                            assert (int(poly[1])*int(f(root)-1)) % gcd(m, int(poly(p))) == 0
                    counts["zero resultants " + str((d1 == 0)+(d2 == 0))] += 1
    return dict(counts)


def example():
    p = 10**45+9
    a, b, c = 112, 292, 312
    powers = [p, p+a, p+b, p+c]
    for w in powers:
        certify_prime(w)
    assert p > 4*c**4+1 and p+b > a-1+(b-1)*b*(b-a)**2*(c-b)
    assert p >= (8*c)**11
    n = prod(powers)
    assert len(str(n)) == len(str(n+1)) == 181


if __name__ == "__main__":
    example()
    print("Certified the four primes and both 181-digit weights.")
    print("Middle-power tests", middle_layer())
    print("Coefficient splits", polynomial_nodes())
    print("Linear bin tests", linear_bins())
