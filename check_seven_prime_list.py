"""Check the optimal four-prime list at 1 + 2*5**4 + 5**52 + 3*5**67."""

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


if __name__ == '__main__':
    check()
