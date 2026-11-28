"""Independent audits for cross-cut ratio bounds. Requires python-flint."""
import argparse
from math import comb, gcd
from pathlib import Path
from flint import fmpz


def factorial_valuation(n, p):
    result = 0
    while n:
        n //= p
        result += n
    return result


def binomial_valuation(n, a, p):
    return (factorial_valuation(n, p) - factorial_valuation(a, p)
            - factorial_valuation(n-a, p))


def proven_prime_factors(n):
    factors = fmpz(n).factor()
    remaining = n
    primes = []
    for p, e in factors:
        assert p.is_prime()
        remaining, rem = divmod(remaining, int(p)**e)
        assert not rem
        primes.append(int(p))
    assert remaining == 1
    return primes


def common_primes(n, a, b):
    for endpoint in (n, n-1):
        for p in proven_prime_factors(endpoint):
            if binomial_valuation(n, a, p) and binomial_valuation(n, b, p):
                yield p


def audit_log(path):
    count = high = largest_witness = 0
    actual, reported = set(), set()
    summary = None
    for line in path.read_text().splitlines():
        if line.startswith("Systems "):
            summary = line
            continue
        tag, *values = line.split()
        x, y, n, a, b = map(int, values)
        assert 2 <= x <= 4096 and x <= y and gcd(x, y) == 1
        assert a*y == b*x and 2 <= a <= b < n
        if tag == "PASS":
            reported.add((n, a, b))
            continue
        assert tag == "ODD"
        assert min(binomial_valuation(n, a, 2), binomial_valuation(n, b, 2)) == 0
        witness = min(common_primes(n, a, b), default=0)
        if n > 10**9:
            assert witness
            high += 1
            largest_witness = max(largest_witness, witness)
        if not witness:
            actual.add((n, a, b))
        count += 1
    assert actual == reported and len(actual) == 162
    assert (count, high, largest_witness) == (99489, 5274, 1427)
    assert summary == ("Systems 7432387, passing parity 99489, passing all bases 162, "
                       "largest weight 36800424470974099")
    print("Directly factored parity survivors:", count, "above one billion:", high,
          "largest least common prime above one billion:", largest_witness)
    return actual


def carry_free_indices(n, p):
    values, place = [0], 1
    while n:
        n, digit = divmod(n, p)
        lower = values
        values = [u+place*d for d in range(digit+1) for u in lower]
        place *= p
    return values


def index_count(n, p):
    count = 1
    while n:
        n, digit = divmod(n, p)
        count *= digit+1
    return count


def row_can_split(n, cut, primes):
    if not primes:
        return True
    for size in (cut, n-cut):
        p = min(primes, key=lambda q: index_count(size, q))
        for i in carry_free_indices(size, p):
            if 1 <= i <= size//2 and all(binomial_valuation(size, i, q) == 0
                                         for q in primes if q != p):
                return True
    return False


def audit_extensions(pairs):
    assignments = 0
    for n, a, b in pairs:
        required_a, required_b, free = [], [], []
        for endpoint in (n, n-1):
            for p in proven_prime_factors(endpoint):
                va = binomial_valuation(n, a, p)
                vb = binomial_valuation(n, b, p)
                assert min(va, vb) == 0
                if va:
                    required_b.append(p)
                elif vb:
                    required_a.append(p)
                else:
                    free.append(p)
        for mask in range(1 << len(free)):
            bases_a = required_a + [p for i,p in enumerate(free) if mask >> i & 1]
            bases_b = required_b + [p for i,p in enumerate(free) if not (mask >> i & 1)]
            assert not (row_can_split(n, a, bases_a) and row_can_split(n, b, bases_b))
            assignments += 1
    return assignments


def divisors(n):
    values = [1]
    for p, e in fmpz(n).factor():
        values = [v * int(p)**j for v in values for j in range(e+1)]
    return values


def small_ratio(limit):
    passed = []
    for x in range(2, limit+1):
        for v in divisors(x):
            for j in range(1, v):
                if gcd(j, v) != 1:
                    continue
                for g in range(1, (4*x*x-1)//27+1):
                    a = g*x
                    product = (a//v)*(a-1)*j*(v-j)
                    for divisor in divisors(product):
                        n = divisor+1
                        if n%v or n <= a:
                            continue
                        b = (n//v)*j
                        if a > b or gcd(a, b) != g:
                            continue
                        if not any(common_primes(n, a, b)):
                            passed.append((x, b//g, n, a, b))
    return sorted(passed)


def exact_pairs(limit):
    count = 0
    for n in range(3, limit+1):
        cuts = [comb(n, a) for a in range(n)]
        for a in range(2, n):
            for b in range(a, n):
                d, j = (gcd(gcd(cuts[a], cuts[b]), v) for v in (n, n-1))
                g = gcd(a, b)
                x, y = a//g, b//g
                da, db = n//gcd(n, a), n//gcd(n, b)
                c = gcd(da, db)
                m, rem = divmod(c*g*x*y, n)
                assert d%c == 0 and not rem and 1 <= m < c*x
                t = c*x-m
                h, rem = divmod(j*m*t*(c*y-m), n-1)
                assert h > 0 and not rem
                k, rem = divmod(g*h, m)
                assert k > 0 and not rem
                ell = j*m*t-k*x
                assert ell > 0 and c*ell*y+h == j*m*m*t
                assert 3125*h*ell*n < 108*j*j*(c*x)**5
                count += 1
    print("Exact binomial pairs:", count)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit-log", type=Path)
    options = parser.parse_args()
    exact_pairs(100)
    pairs = small_ratio(32)
    assert pairs == [(10, 23, 46, 10, 23), (23, 36, 46, 23, 36),
                     (23, 43, 989, 92, 172)]
    print("Endpoint-coprime pairs through reduced term 32:", pairs)
    if options.audit_log:
        survivors = audit_log(options.audit_log)
        assert len(survivors) == 162 and max(survivors)[0] == 55784447
        assignments = audit_extensions(survivors)
        assert assignments == 164
        print("Unextendible cut pairs:", len(survivors), "endpoint-prime assignments:", assignments)
