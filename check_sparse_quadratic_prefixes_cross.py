"""Check the finite quadratic range by weight prefixes alone and intervals for P."""
from itertools import product
from math import isqrt

from check_sparse_quadratic_powers import KNOWN_SCALE, RESIDUAL_ROWS


def allowed_prefix(a, slopes, d):
    # The discrepancy inequality is equivalent to
    # s*(a+1)/d - 1 <= v <= s*(a+1)/d in each coordinate.
    floors = []
    flexible = []
    for slot, s in enumerate(slopes):
        v, rem = divmod(s*(a+1), d)
        floors.append(v)
        if rem == 0 and v:
            flexible.append(slot)
    for corrections in product((0, 1), repeat=len(flexible)):
        v = floors.copy()
        for slot, change in zip(flexible, corrections):
            v[slot] -= change
        if sum(v) == a and not (v[0] & v[1] or v[0] & v[2] or v[1] & v[2]):
            return True
    return False


def check_shape(row, limit):
    d, i, j = row
    slopes = (3*i, 3*j, d-3*(i+j))
    frontier = {0}
    count = 0
    shift = (d-1).bit_length()-1
    for h in range(1, limit.bit_length()-shift):
        frontier = {v for a in frontier for v in (2*a, 2*a+1)
                    if v and allowed_prefix(v, slopes, d)}
        b = h+shift
        lo = max(1 << b, KNOWN_SCALE+1)
        hi = min((2 << b)-1, limit)
        if lo > hi:
            continue
        for n in range((lo*(3*lo+2)).bit_length(), (hi*(3*hi+2)).bit_length()+1):
            for a in frontier:
                left = 3*a << (n-h)
                right = 3*(a+1) << (n-h)
                start = max(lo, (isqrt(left)+2)//3)
                stop = min(hi, (isqrt(right-2)-1)//3)
                first = start+(pow(3, -1, d)-start) % d
                for p in range(first, stop+1, d):
                    q = (3*p-1)//d
                    count += 1
                    if q % 2:
                        v = 3*p+2
                        parts = (v*q*i, v*q*j, v*(p-q*(i+j)))
                        assert sum(bin(x).count('1') for x in parts) > bin(p*v).count('1')
    return count


if __name__ == "__main__":
    count = sum(check_shape(row, 10**1000) for row in RESIDUAL_ROWS)
    assert count == 120253
    print("Integral inputs examined by the scalar prefix check:", count)
