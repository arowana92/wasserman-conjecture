"""Check optimal four-prime lists at two seven-unit weights."""

from itertools import product

from check_dense_offsets import free_rows
from check_scaled_row_blocks import certify_prime
from check_three_factors import valuation


def check():
    n = 1 + 2*5**4 + 5**52 + 3*5**67
    q, s, r = n//31, (n-1)//(5**4*5922), n-38
    for prime in (5, 31, q, s, r):
        certify_prime(prime)
    assert n == 31*q and q != 31 and n-1 == 5**4*5922*s
    assert len({5, q, s, r}) == 4 and n % s == 1
    assert q*(2*s+1) > n and r > n//2 and n-r < min(s, 2*q)
    rows = set(free_rows(n, 5))
    # Independently distribute the seven labelled powers instead of splitting digits.
    units = (1, 5**4, 5**4, 5**52, 5**67, 5**67, 5**67)
    labelled = set()
    for labels in product(range(3), repeat=len(units)):
        parts = [0, 0, 0]
        for power, label in zip(units, labels):
            parts[label] += power
        if all(parts):
            labelled.add(tuple(sorted(parts)))
    assert rows == labelled and len(rows) == 68
    for row in rows:
        assert sum(row) == n and valuation(n, row, 5) == 0
        assert all(valuation(n, row, prime) > 0 for prime in (q, s, r))
    # A direct independent enumeration of every possible q-free row.
    for a in range(1, 31//3+1):
        for b in range(a, (31-a)//2+1):
            row = (q*a, q*b, q*(31-a-b))
            assert all(valuation(n, row, prime) > 0 for prime in (5, s, r))
    assert n % 2 == 1 and (n-1) % 5 == 0
    print(f'N = {n}\nq = {q}\ns = {s}\nr = {r}')
    print('68 carry-free rows at 5. The four-prime certificate is optimal.')


def check_distinct_levels():
    p = 257
    places = (1, 2, 3, 7, 17, 23)
    n = 1 + sum(p**h for h in places)
    q, s, r = n//3, (n-1)//(1598*p), n-44
    for prime in (3, p, q, s, r):
        certify_prime(prime)
    assert n == 3*q and q != 3 and n-1 == 1598*p*s and len({p, q, s, r}) == 4
    assert n-p**23 < min(s, 2*q) and n-r < 1+p
    assert r > n//2 and n-r < min(2*q, s) and q*(2*s+1) > n
    rows = set(free_rows(n, p))
    assert len(rows) == 301 and n % 2 == 1 and (n-1) % p == 0
    for row in rows:
        assert valuation(n, row, p) == 0
        assert all(valuation(n, row, t) > 0 for t in (q, s, r))
    print(f'Four-prime certificate at N = {n} with {len(rows)} carry-free rows at 257.')


if __name__ == '__main__':
    check()
    check_distinct_levels()
