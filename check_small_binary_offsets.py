"""Check all remaining exponents in the binary shape theorem at offsets through 256."""

from collections import Counter

from check_dense_offsets import free_rows, has_size_obstruction


def low_pairs(d):
    # With a zero difference the lower contributions are b, 2b and d-3b.
    return [(b, d-3*b) for b in range(1, d//3+1)
            if (b | (2*b) | (d-3*b)) == d]


def audit():
    counts = Counter()
    for d in range(1, 257):
        pairs = low_pairs(d)
        for k in range(2, (d**25-1).bit_length()):
            q = 1 << k
            if 4*q+d <= 10**9:
                continue
            for ell in range(1, k):
                middle = 1 << ell
                n = 3*q+middle+d
                if n <= 10**9:
                    continue
                counts['weight and exponent choices'] += 1
                threshold = 27*(n-1)**2
                assert (2*d)**5 < threshold and 2*d < q
                if middle <= d:
                    continue
                counts['three separate high bits'] += 1
                assert (2*d+2)**7 < threshold
                r = q//middle
                for b, f in pairs:
                    counts['zero differences'] += 1
                    if b == r*f:
                        counts['permitted unshifted rows'] += 1
                    else:
                        row = (q+b, 2*(q+b), middle+f)
                        assert has_size_obstruction(n, row), (n, d, row)
    print('Remaining exponent scan:', dict(counts))


def sample_control():
    cases = rows = passing = 0
    for k in (31, 53, 79):
        for ell in (2, k//2, k-3, k-1):
            q, middle, r = 1 << k, 1 << ell, 1 << (k-ell)
            for d in range(1, 65):
                if d**25 <= q:
                    continue
                n = 3*q+middle+d
                cases += 1
                seen = set()
                for row in free_rows(n, 2):
                    rows += 1
                    if not has_size_obstruction(n, row):
                        seen.add(row)
                expected = set()
                u, rem = divmod(d, 3*r+1)
                if rem == 0 and (u | r*u | 2*r*u) == d:
                    f = middle+u
                    row = (f, r*f, 2*r*f)
                    if not has_size_obstruction(n, row):
                        expected.add(row)
                assert seen == expected
                passing += len(seen)
    assert passing
    print('Direct control weights, rows and passing rows:', cases, rows, passing)


if __name__ == '__main__':
    audit()
    sample_control()
