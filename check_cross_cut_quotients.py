"""Enumerate pairs by lcm(a,b)/N using a symmetric divisor equation."""
import argparse
from math import gcd, isqrt, lcm
from pathlib import Path
from check_cross_cut_audit import audit_extensions, common_primes, divisors


def small_factors(value, step, residue, low):
    high = isqrt(value)
    first = low+(residue-low)%step
    if first > high:
        return []
    if (high-first)//step < 45:
        return (d for d in range(first,high+1,step) if value%d == 0)
    return (d for d in divisors(value) if low <= d <= high and (d-residue)%step == 0)


def scan(limit, cutoff):
    systems = parity = largest = 0
    passed = set()
    for m in range(2, limit+1):
        square = m*m
        # The smaller of x,y is below 2*m^2/u.
        last = min(m-1, (2*square-1)//(cutoff+1))
        for u in range(1,last+1):
            k = m-u
            low = max(1, u*(cutoff+1)-square)
            for h in divisors(k*m):
                value = k*m*square-u*h
                if value < low*low:
                    continue
                for ell in small_factors(value, u, -square, low):
                    x = (square+ell)//u
                    y, rem = divmod(square+value//ell, u)
                    if rem or x <= m or gcd(x,y) != 1 or h%gcd(m,x) or k*x*y%h:
                        continue
                    n, g = k*x*y//h, k*m//h
                    a, b = g*x, g*y
                    assert cutoff < x <= y and 2 <= a <= b < n and lcm(a,b) == m*n
                    systems += 1
                    largest = max(n, largest)
                    if a & (n-a) and b & (n-b):
                        continue
                    parity += 1
                    if not any(common_primes(n,a,b)):
                        passed.add((n,a,b))
                        print("PASS", m,x,y,n,a,b, flush=True)
        if m%100 == 0:
            print("Through quotient", m, "systems", systems, "parity", parity, flush=True)
    return systems, parity, largest, passed


def audit_cpp(path):
    systems = parity = largest = None
    passed, seen = set(), set()
    for line in path.read_text().splitlines():
        if line.startswith("Systems "):
            _, s, _, p, _, n = line.split()
            systems, parity, largest = int(s), int(p), int(n)
            continue
        n, a, b = map(int, line.split())
        assert 2 <= a <= b < n and lcm(a,b)%n == 0
        assert not (a & (n-a)) or not (b & (n-b))
        record = (n,a,b)
        assert record not in seen
        seen.add(record)
        if not any(common_primes(n,a,b)):
            passed.add(record)
    assert len(seen) == parity
    return systems, parity, largest, passed


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("limit", type=int)
    parser.add_argument("cutoff", type=int, help="require the smaller ratio term to exceed this")
    parser.add_argument("--audit-log", type=Path)
    options = parser.parse_args()
    systems, parity, largest, passed = scan(options.limit, options.cutoff)
    print("Systems:", systems, "parity:", parity, "largest weight:", largest)
    assignments = audit_extensions(passed)
    print("Unextendible pairs:", len(passed), "assignments:", assignments)
    if (options.limit, options.cutoff) == (1000, 4096):
        assert (systems, parity, largest, len(passed), assignments) == (
            8102110, 38222, 998002997000999999001, 19, 20)
    if options.audit_log:
        assert audit_cpp(options.audit_log) == (systems, parity, largest, passed)
        print("Independent enumeration agrees.")
