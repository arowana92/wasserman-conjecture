"""Check the finite prime-power bases by direct subset sums and factorial floors."""

from argparse import ArgumentParser
from collections import defaultdict
from itertools import combinations, product
from math import comb, gcd, isqrt, prod

from flint import fmpz


def prime_powers(limit):
    result = set()
    for p in range(2, limit + 1):
        if any(p % t == 0 for t in range(2, isqrt(p) + 1)):
            continue
        q = p
        while q <= limit:
            result.add(q)
            q *= p
    return sorted(result)


def factor(value):
    factors = [(p, int(e)) for p, e in fmpz(value).factor()]
    assert prod(p**e for p, e in factors) == value
    assert all(p.is_prime() for p, _ in factors)
    return [(int(p), int(p**e)) for p, e in factors]


def coefficients_avoid(row, prime):
    power = prime
    total = sum(row)
    while power <= total:
        if total // power != sum(a // power for a in row):
            return False
        power *= prime
    return True


def exclude_pair(weight, factors):
    """Enumerate the row receiving S >= V by its three exceptional positions."""
    predecessor = factor(weight - 1)
    for take_weight in product((False, True), repeat=len(factors)):
        upper = prod(w for (_, w), take in zip(factors, take_weight) if take)
        lower = weight // upper
        if min(upper, lower) < 3:
            continue
        for take_shift in product((False, True), repeat=len(predecessor)):
            left = [item for item, take in zip(predecessor, take_shift) if take]
            s = prod(w for _, w in left)
            if s < lower or 27 * ((weight - 1) // s) ** 2 > upper**3:
                continue
            required = [p for (p, _), take in zip(factors, take_weight) if take]
            required.extend(p for p, _ in left)
            for positions in product(range(3), repeat=len(left)):
                if len(set(positions)) < 3:
                    continue
                parts = []
                for i in range(3):
                    zero = prod(w for (_, w), j in zip(left, positions) if j != i)
                    one = s // zero
                    parts.append(upper * zero * (lower * pow(zero, -1, one) % one))
                if sum(parts) == weight and all(coefficients_avoid(parts, p) for p in required):
                    raise ArithmeticError(("heavy row survives", weight, parts))
    return True


def subsets(powers, max_size):
    return [[sum(part) for part in combinations(powers, k)] for k in range(max_size + 1)]


def necklaces(n, limit):
    """Least binary cyclic rotations, in lexicographic order."""
    digits = [0] * (n + 1)

    def extend(i, period, used, word):
        if i > n:
            if n % period == 0 and used >= 2:
                yield word
            return
        first = digits[i - period]
        if used + first <= limit:
            digits[i] = first
            yield from extend(i + 1, period, used + first, 2 * word + first)
        if first == 0 and used < limit:
            digits[i] = 1
            yield from extend(i + 1, i, used + 1, 2 * word + 1)

    yield from extend(1, 1, 0, 0)


def prime_catalog(n):
    k = n // 3
    # The norm is positive by conjugate pairing. One splitting prime larger
    # than its Parseval bound recovers its exact integer value.
    limit = isqrt((k * (n - k)) ** (n - 1) // (n - 1) ** (n - 1)) + 1
    modulus = 1 + n * (limit // n + 1)
    while not fmpz(modulus).is_prime():
        modulus += n
    z = 1
    seed = 1
    while z == 1:
        seed += 1
        z = pow(seed, (modulus - 1) // n, modulus)
    assert pow(z, n, modulus) == 1
    values = [[pow(z, j * i, modulus) for i in range(n)] for j in range(1, n)]
    full = (1 << n) - 1
    seen = set()
    primes = set()
    count = 0
    for word in necklaces(n, k):
        count += 1
        if word in seen:
            continue
        positions = [i for i in range(n) if word & (1 << i)]
        for h in range(1, n):
            image = sum(1 << ((h * i) % n) for i in positions)
            least = image
            for _ in range(n - 1):
                image = ((image << 1) & full) | (image >> (n - 1))
                least = min(least, image)
            seen.add(least)
        residue = 1
        for row in values:
            residue = residue * sum(row[i] for i in positions) % modulus
        assert 1 <= residue < limit
        primes.update(s for s, _ in factor(residue) if (s - 1) % n == 0)
    assert len(seen) == count
    assert n * count == sum(comb(n, i) for i in range(2, k + 1))
    return prod(primes)


def forced_power(n, q, upper):
    """Find a predecessor power too large for the complementary normalized row."""
    from flint import fmpz_poly

    predecessor = (q**n - 1) // (q - 1) - 1
    best = 1
    for d in range(2, n):
        if (n - 1) % d:
            continue
        poly = fmpz_poly.cyclotomic(d)
        if poly.degree() > 6:
            continue
        for p, w in factor(int(poly(q))):
            while predecessor % (w * p) == 0:
                w *= p
            if 2 * w + 1 > upper:
                best = max(best, w)
    return best


def verify(n, q, catalog):
    weight = (q**n - 1) // (q - 1)
    factors = factor(gcd(weight, catalog)) if catalog is not None else factor(weight)
    max_size = n // 3
    # A prime dividing a part and R_n(q) divides its algebraic norm.
    # Parseval's upper bound increases with k up to n/3.
    possible = [
        p
        for p, _ in factors
        if p**2 * (n - 1) ** (n - 1) <= (max_size * (n - max_size)) ** (n - 1)
    ]
    if not possible:
        return False
    upper = 1
    for p in possible:
        remaining = weight
        while remaining % p == 0:
            upper *= p
            remaining //= p
    w = forced_power(n, q, upper) if catalog is not None else 1
    powers = [q**i for i in range(n)]
    a = subsets(powers[: n // 2], max_size)
    b = subsets(powers[n // 2 :], max_size)
    for p in possible:
        modulus = p * w
        targets = {0, p * pow(p, -1, w) % modulus} if w > 1 else {0}
        for j, parts in enumerate(b):
            index = defaultdict(list)
            for value in parts:
                for target in targets:
                    index[(target - value) % modulus].append(value)
            for k in range(max(0, 2 - j), max_size - j + 1):
                for x in a[k]:
                    for y in index.get(x % modulus, ()):
                        f = x + y
                        u = gcd(weight, f)
                        t = (weight - 1) // gcd(weight - 1, f * (f - 1))
                        if 27 * t * t <= u**3:
                            return exclude_pair(weight, factor(weight))
    return False


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("length", type=int)
    parser.add_argument("cutoff", type=int)
    parser.add_argument(
        "--catalog",
        action="store_true",
        help="calculate all possible weight primes by finite-field norms",
    )
    args = parser.parse_args()
    assert args.length >= 3 and all(args.length % i for i in range(2, isqrt(args.length) + 1))
    bases = prime_powers(args.cutoff)
    catalog = prime_catalog(args.length) if args.catalog else None
    fallback = [q for q in bases if verify(args.length, q, catalog)]
    print(args.length, args.cutoff, len(bases), "bases checked, CRT needed at", fallback)
