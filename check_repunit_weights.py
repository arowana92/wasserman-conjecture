"""Certify prime-length repunit weights at every prime-power base."""

from argparse import ArgumentParser
from functools import lru_cache
from itertools import combinations
from math import comb, gcd, isqrt, prod

from flint import fmpz, fmpz_mat
from search_endpoint import find_row, splits


@lru_cache(None)
def prime(r):
    return r >= 2 and all(r % i for i in range(2, isqrt(r) + 1))


def floor_root(x, n):
    lo, hi = 0, 1 << ((x.bit_length() + n - 1) // n)
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if mid**n <= x:
            lo = mid
        else:
            hi = mid
    return lo


def orbits(n, k):
    """Normalize two ordered positions to 0 and 1 by an affine map mod n."""
    seen = set()
    for tail in combinations(range(2, n), k - 2):
        points = (0, 1, *tail)
        mask = sum(1 << j for j in points)
        if mask in seen:
            continue
        stabilizer = 0
        for first in points:
            for second in points:
                if first == second:
                    continue
                inv = pow(second - first, -1, n)
                image = sum(1 << ((j - first) * inv % n) for j in points)
                seen.add(image)
                stabilizer += image == mask
        yield points, n * (n - 1) // stabilizer


def norm(n, points):
    # Multiplication by F in the basis 1,X,...,X^(n-2) of Z[X]/R_n.
    # First use X^n=1, then X^(n-1)=-1-X-...-X^(n-2).
    bits = [int(i in points) for i in range(n)]
    return abs(
        int(
            fmpz_mat(
                [
                    [bits[(i - j) % n] - bits[(n - 1 - j) % n] for j in range(n - 1)]
                    for i in range(n - 1)
                ]
            ).det()
        )
    )


@lru_cache(None)
def count_bits(word):
    return bin(word).count("1")


def can_split_complement(n, points, root, modulus):
    """Can a second part of at least k units vanish at the same root?"""
    k = len(points)
    remaining = [pow(root, i, modulus) for i in range(n) if i not in points]
    mid = len(remaining) // 2
    left = half_sums(remaining[:mid], modulus)
    right = half_sums(remaining[mid:], modulus)
    sizes = {}
    for mask, value in enumerate(left):
        sizes[value] = sizes.get(value, 0) | (1 << count_bits(mask))
    for mask, value in enumerate(right):
        lo = max(0, k - count_bits(mask))
        hi = n - 2 * k - count_bits(mask)
        if hi >= lo and (sizes.get(-value % modulus, 0) >> lo) & ((1 << (hi - lo + 1)) - 1):
            return True
    return False


def weight_primes(n, points):
    absolute = norm(n, points)
    assert absolute > 0
    factors = [(int(r), int(v)) for r, v in fmpz(absolute).factor()]
    assert prod(r**v for r, v in factors) == absolute and all(prime(r) for r, _ in factors)
    support = []
    for r, v in factors:
        if (r - 1) % n:
            continue
        g = 2
        while pow(g, (r - 1) // n, r) == 1:
            g += 1
        z = pow(g, (r - 1) // n, r)
        allowed = set()
        power = 1
        for j in range(1, v + 1):
            w = r**j
            root = pow(z, r ** (j - 1), w)
            assert root % r != 1 and pow(root, n, w) == 1
            roots = {
                pow(root, h, w)
                for h in range(1, n)
                if sum(pow(root, h * i, w) for i in points) % w == 0
            }
            if w > 1000:
                roots = {z for z in roots if can_split_complement(n, points, z, w)}
            if not roots:
                break
            power = w
            if j == 1:
                allowed = roots
        if allowed:
            support.append((r, power, allowed))
    return support


def half_sums(powers, modulus):
    sums = [0]
    for value in powers:
        sums += [(s + value) % modulus for s in sums]
    return sums


class BadRoots:
    """Lower-bound the number of obstructing complex roots by exact modular tests."""

    def __init__(self, n):
        m = n - 1
        modulus = 1 + 1000 * m
        while not prime(modulus):
            modulus += m
        self.half = n // 2
        self.low_mask = (1 << self.half) - 1
        self.tables = []
        for d in range(2, m + 1):
            if m % d:
                continue
            units = [a for a in range(1, d) if gcd(a, d) == 1]
            g = 2
            while True:
                w = pow(g, (modulus - 1) // d, modulus)
                if all(pow(w, j, modulus) != 1 for j in range(1, d)):
                    break
                g += 1
            powers = [pow(w, i, modulus) for i in range(n)]
            self.tables.append(
                (
                    len(units),
                    half_sums(powers[: self.half], modulus),
                    half_sums(powers[self.half :], modulus),
                )
            )
        self.modulus = modulus

    def degree(self, word):
        lo, hi = word & self.low_mask, word >> self.half
        return sum(
            degree for degree, a, b in self.tables if (a[lo] + b[hi]) % self.modulus not in (0, 1)
        )


@lru_cache(maxsize=4096)
def powers_at(n, q):
    powers = tuple(q**j for j in range(n))
    return sum(powers), powers


def cut_passes(n, q, mask):
    weight, powers = powers_at(n, q)
    f = 0
    if q == 2:
        f = mask
    else:
        while mask:
            bit = mask & -mask
            f += powers[bit.bit_length() - 1]
            mask -= bit
    u = gcd(f, weight)
    t = (weight - 1) // gcd(weight - 1, f * (f - 1))
    return u >= 3 and 27 * t * t <= u**3


def full_powers(value):
    factors = fmpz(value).factor()
    assert all(p.is_prime() for p, _ in factors)
    powers = [(int(p), int(p**v)) for p, v in factors]
    assert prod(w for _, w in powers) == value
    return powers


@lru_cache(None)
def exclude_weight(n, q):
    weight, _ = powers_at(n, q)
    weight_factors = full_powers(weight)
    predecessor_splits = list(splits(full_powers(weight - 1)))
    for u, first, second in splits(weight_factors):
        v = weight // u
        if min(u, v) < 3:
            continue
        for s, left, right in predecessor_splits:
            t = (weight - 1) // s
            # At least one of S >= V and T >= U holds, by integrality.
            if s < v or len(left) < 3 or 27 * t * t > u**3:
                continue
            if find_row(weight, u, first, s, left) and find_row(weight, v, second, t, right):
                raise ArithmeticError(("unresolved pair of rows", n, q, u, s))
    return True


def check_size(n, k, bad):
    masks_checked = candidates = visited = largest = largest_prime = 0
    full = (1 << n) - 1
    m = n - 1
    for points, orbit_size in orbits(n, k):
        visited += orbit_size
        support = weight_primes(n, points)
        if not support:
            continue
        largest_prime = max(largest_prime, *(r for r, _, _ in support))
        upper_u = prod(power for _, power, _ in support)
        upper_t = isqrt(upper_u**3 // 27)
        seen = set()
        for h in range(1, n):
            base = sum(1 << ((h * i) % n) for i in points)
            inverse = pow(h, -1, n)
            residues = [(r, {pow(s, inverse, r) for s in roots}) for r, _, roots in support]
            base_cache = {}
            for shift in range(n):
                word = ((base << shift) | (base >> (n - shift))) & full
                if word in seen:
                    continue
                seen.add(word)
                e = bad.degree(word)
                # Cyclic coefficient comparison is also valid when the modular
                # calculation certifies no roots for this particular word.
                limit = 1 + (k + 4) * upper_t
                if e:
                    c = min(k * (k + 1), (m * (k + 3) + e - 1) // e)
                    limit = min(limit, 1 + floor_root(upper_t * c**e, e))
                largest = max(largest, limit)
                if limit not in base_cache:
                    base_cache[limit] = {
                        q
                        for r, roots in residues
                        for s in roots
                        for q in range(s, limit + 1, r)
                        if q >= 2
                    }
                for q in base_cache[limit]:
                    candidates += 1
                    if cut_passes(n, q, word):
                        exclude_weight(n, q)
        assert len(seen) == orbit_size
        masks_checked += len(seen)
    assert visited == comb(n, k)
    return visited, masks_checked, largest, candidates, largest_prime


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument(
        "lengths", nargs="*", type=int, default=[3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
    )
    parser.add_argument("--size", type=int, help="only one size of the smallest digit subset")
    args = parser.parse_args()
    for n in args.lengths:
        if n < 3 or not prime(n):
            parser.error("each length must be an odd prime")
        if args.size is not None and not 2 <= args.size <= n // 3:
            parser.error("subset size must be between 2 and length//3")
        bad = BadRoots(n)
        sizes = [args.size] if args.size is not None else range(2, n // 3 + 1)
        results = [check_size(n, k, bad) for k in sizes]
        print(
            n,
            sum(r[0] for r in results),
            sum(r[1] for r in results),
            max((r[2] for r in results), default=0),
            sum(r[3] for r in results),
            max((r[4] for r in results), default=0),
            exclude_weight.cache_info().misses,
            flush=True,
        )
        exclude_weight.cache_clear()
