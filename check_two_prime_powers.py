"""Check the four-bin cofactor condition and its implications for individual prime powers."""

from itertools import product
from math import comb, gcd, isqrt, lcm, prod

from check_endpoint_pair_bound import prime_powers, rows


def check_partitions():
    for d in range(1, 5001):
        factors = list(prime_powers(d))
        possible = False
        for colors in product(range(4), repeat=len(factors)):
            bins = [1] * 4
            for (_, power), place in zip(factors, colors):
                bins[place] *= power
            a, b, c, e = bins
            if (a+b)*(c+e) < d:
                possible = True
                break
        excluded = (len(factors) <= 1 or d in (12, 15)
                    or (len(factors) == 2 and factors[0] == (2, 2)))
        assert possible != excluded, d
    print("Four-bin classification checked through 5000.")


def avoiding_rows(n, scale, factors):
    table = []
    for normalized, _ in rows(n//scale):
        x, y, _ = normalized
        low = {w for _, w in factors if sum(v % w == 0 for v in normalized) == 2}
        if not low:
            continue
        coefficient = comb(n, scale*x)*comb(n-scale*x, scale*y)
        avoided = {w for p, w in factors if coefficient % p}
        assert avoided <= low
        for w in low:
            zeros = [v//w for v in normalized if v % w == 0]
            table.append((w, zeros, low, avoided))
    return table


def check_assignments():
    cases = observed = actual = 0
    for n in range(3, 1801):
        factors = list(prime_powers(n-1))
        if len(factors) < 2:
            continue
        for u in range(3, isqrt(n)+1):
            if n % u or gcd(u, n//u) != 1:
                continue
            v = n//u
            if v > 100:
                continue
            left = avoiding_rows(n, u, factors)
            right = avoiding_rows(n, v, factors)
            for q, (x, y), allowed_a, exact_a in left:
                for r, (z, t), allowed_c, exact_c in right:
                    if q == r:
                        continue
                    d = (n-1)//(q*r)
                    assert (x+y)*(z+t) < d
                    covered = (allowed_a | allowed_c) - {q, r}
                    assert lcm(x, y, z, t) % prod(covered) == 0
                    cases += 1
                    observed += bool(covered)
                    actual += q in exact_a and r in exact_c
    assert cases and observed and actual
    print(f"Checked {cases} pairs of low-residue assignments, {observed} covering a further power, {actual} actually avoiding both assigned primes.")


if __name__ == "__main__":
    check_partitions()
    check_assignments()
