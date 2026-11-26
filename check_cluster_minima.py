"""Certify the initial interval of optimal four-prime lists at a 139-digit product."""

from concurrent.futures import ProcessPoolExecutor
from math import gcd, prod

from flint import fmpz

from check_cluster_certificates import is_prime, universal_prime


CORE = tuple(10**46 + offset for offset in (121, 537, 117073))
PRODUCT = prod(CORE)
LAST = 1123


def integer_root(n, exponent):
    """Floor root, using only integer comparisons."""
    lo = 1 << ((n.bit_length()-1)//exponent)
    hi = 2*lo
    while hi-lo > 1:
        mid = (hi+lo)//2
        if mid**exponent <= n:
            lo = mid
        else:
            hi = mid
    return lo


def candidate_bases(n):
    # Include composite bases so the digit check never requires factoring n.
    yield n
    yield n-1
    if n % 2 == 0:
        yield n//2
    for j in range(2, n.bit_length()):
        b = integer_root(n, j)
        left = n-b**j
        if left == 0:
            yield b
        else:
            while left % b == 0:
                left //= b
            if left == 1:
                yield b
        if n % 2 == 0 and 2**j <= n//2:
            b = integer_root(n//2, j)
            if 2*b**j == n:
                yield b


def composite_witness(n):
    for a in range(2, 100):
        if a >= n:
            return False
        if gcd(a, n) > 1 or pow(a, n-1, n) != 1:
            return True
    return False


def check_pair(d):
    n = d*PRODUCT
    for weight in (n, n+1):
        assert universal_prime(weight) == 0
        for b in set(candidate_bases(weight)):
            assert composite_witness(b), (d, weight, b)
    fourth = n-1 if n % 2 == 0 else n-2
    while not (fmpz(fourth).is_probable_prime() and is_prime(fourth)):
        fourth -= 2
    assert n+1-fourth <= CORE[0] and 2*fourth > n+1
    return n-fourth, d


def main():
    p, q, r = CORE
    assert all(is_prime(a) for a in CORE)
    assert r-q > q-p and 2*(LAST+1)*(r-p)**2 <= p
    # Small inputs validate the exhaustive list of bases against digit division.
    controls = 0
    for n in range(4, 350):
        found = set(candidate_bases(n))
        for b in range(2, n+1):
            k, digit_sum = n, 0
            while k:
                k, digit = divmod(k, b)
                digit_sum += digit
            assert (b in found) == (digit_sum <= 2)
            controls += 1
    with ProcessPoolExecutor(max_workers=8) as pool:
        largest_gap = max(pool.map(check_pair, range(1, LAST+1), chunksize=8))
    boundary = (LAST+1)*PRODUCT
    assert is_prime(boundary-1) and universal_prime(boundary) == boundary-1
    print(f'{controls} exact base checks; both minimum sizes are four for d=1..{LAST}.')
    print(f'Largest gap to the fourth prime and its multiplier: {largest_gap}.')
    print(f'At d={LAST+1} the weight dM has universal prime dM-1.')


if __name__ == '__main__':
    main()
