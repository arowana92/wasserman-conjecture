"""Reconstruct all assignments with unit multipliers at a bounded companion quotient."""

from argparse import ArgumentParser
from collections import Counter
from functools import lru_cache
from itertools import combinations
from math import gcd, prod

from check_scaled_row_blocks import certify_prime

from flint import fmpz

from check_three_factors import valuation


@lru_cache(maxsize=8192)
def factors(n):
    result = tuple((int(p), int(e)) for p, e in fmpz(n).factor())
    assert all(fmpz(p).is_prime() for p, _ in result) and prod(p**e for p, e in result) == n
    return result


def divisors(n):
    divs = [1]
    for p, e in factors(n):
        divs = [v*p**i for v in divs for i in range(e+1)]
    return sorted(divs)


def companion_bins(limit):
    for j in range(1, limit+1):
        # d_i odd and T=1+jU force j even and U odd.
        if j%2:
            continue
        for r in range(j+1, 3*j+2):
            if r%2 == 0:
                continue
            delta = r-j
            value = (j*r)**2+delta
            for v in divisors(value):
                if v*v > value:
                    break
                s, sr = divmod(v+j*r, delta)
                t, tr = divmod(value//v+j*r, delta)
                if sr or tr or not r <= s <= t or s%2 == 0 or t%2 == 0:
                    continue
                u = r*s+r*t+s*t
                assert r*s*t == 1+j*u and gcd(r*s,r*t,s*t) == 1
                yield j, (r,s,t), u


def compare_small(limit):
    cap = 110
    by_divisors = {(j, bins, u) for j, bins, u in companion_bins(limit) if bins[2] <= cap}
    by_bins = set()
    for r, s, t in combinations(range(3, cap+1, 2), 3):
        total, m = r*s+r*t+s*t, r*s*t
        j, rem = divmod(m-1, total)
        if not rem and 1 <= j <= limit:
            by_bins.add((j, (r,s,t), total))
    assert by_divisors == by_bins and by_bins
    return len(by_bins)


def compare_small_quotients(limit):
    # Solve the original linear equation for t with a direct (r,s) loop.
    # Ordering s <= t gives (r-j)*s <= 2*j*r+1.
    limit = min(limit, 50)
    result = set()
    for j in range(1, limit+1):
        for r in range(j+1, 3*j+2):
            for s in range(r, (2*j*r+1)//(r-j)+1):
                denominator = (r-j)*s-j*r
                if denominator <= 0:
                    continue
                t, remainder = divmod(j*r*s+1, denominator)
                if remainder or s > t or any(k%2 == 0 for k in (r,s,t)):
                    continue
                result.add((j, (r,s,t), r*s+r*t+s*t))
    assert set(companion_bins(limit)) == result
    return len(result)


def example():
    primes = (1000000007, 2147483647, 4294967291, 154742504079676789)
    for p in primes:
        certify_prime(p)
    n = prod(primes)+1
    q = 2**150
    h = n-q
    assert n == 1427247692705961760186075461832537708033263872
    assert q < n < 2*q and bin(h).count('1') == 50
    assert q*q < h**3 and h*h*(h-1) <= 16*(n-1)**2
    assert h**7 > 6**6*(n-1)**4
    print('Certified four predecessor primes at', n)


def check(limit):
    counts = Counter()
    largest_candidate = 0
    for j, bins, u in companion_bins(limit):
        t = prod(bins)
        counts['companion triples'] += 1
        assert u%8 == 3 and (sum(bins)-9*j)%8 == 3
        assert min(bins) > 1
        if all(len(factors(w)) == 1 for w in bins):
            counts['three prime powers'] += 1
            continue
        # (d-U)(e-U)=U^2-1 and one of d,e is even.
        for f in divisors(u*u-1):
            if f*f > u*u-1:
                break
            d, e = u+f, u+(u*u-1)//f
            if (d-e)%2 == 0 or gcd(d,e) != 1 or gcd(d*e, t) != 1:
                continue
            n = d*e*t+1
            v, rem = divmod(n, u)
            assert rem == 0 and u*(d+e) == d*e+1
            if gcd(u, v) != 1:
                continue
            counts['unitary assignments'] += 1
            largest_candidate = max(largest_candidate, n)
            row = (u*d, u*e, u*j*d*e)
            companion = tuple(v*(t//w) for w in bins)
            assert sum(row) == sum(companion) == n
            if row[0]&row[1]:
                continue
            counts['first pair binary carry-free'] += 1
            if (row[0]|row[1])&row[2]:
                continue
            counts['binary carry-free assignments'] += 1
            if any(valuation(n, row, p) for p, _ in factors(u)):
                continue
            counts['also avoids its weight primes'] += 1
            if any(valuation(n, row, p) for k in (d,e) for p, _ in factors(k)):
                continue
            counts['first row takes all its primes'] += 1
            if not any(valuation(n, companion, p) for k in (*bins, v) for p, _ in factors(k)):
                raise ArithmeticError((n, row, companion))
    print('Companion comparison through largest bin 110:', compare_small(limit))
    print('Direct comparison through quotient', min(50, limit), ':', compare_small_quotients(limit))
    print('Upper quotient', limit, dict(counts))
    print('Largest unitary assignment tested:', largest_candidate)


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('limit', type=int, nargs='?', default=1000)
    example()
    check(parser.parse_args().limit)
