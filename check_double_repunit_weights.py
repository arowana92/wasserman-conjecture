"""Certify repunit weights of twice-prime length by cyclotomic norms and CRT."""

from argparse import ArgumentParser
from functools import lru_cache
from itertools import combinations, product
from math import comb, gcd, lcm, prod

from flint import fmpz, fmpz_poly
from check_repunit_weights import floor_root, half_sums, powers_at


@lru_cache(None)
def factorization(value):
    factors = tuple((int(p), int(e)) for p, e in fmpz(value).factor())
    assert prod(p**e for p, e in factors) == value
    assert all(fmpz(p).is_prime() for p, _ in factors)
    return factors


def divisors(n):
    return [d for d in range(2, n + 1) if n % d == 0]


@lru_cache(None)
def roots_of_order(order, prime, exponent):
    """Roots of the cyclotomic polynomial modulo a candidate prime power."""
    h = order
    while h % prime == 0:
        h //= prime
    if h != order and exponent > 1:
        # Only odd orders occur in this local refinement. At an odd prime
        # dividing the order, a cyclotomic value has valuation at most one.
        assert prime != 2
        return ()
    if (prime - 1) % h:
        return ()
    if h == 1:
        return (1,)
    factors = [int(p) for p, _ in factorization(h)]
    seed = 2
    while True:
        z = pow(seed, (prime - 1) // h, prime)
        if all(pow(z, h // p, prime) != 1 for p in factors):
            break
        seed += 1
    modulus = prime**exponent
    root = pow(z, prime ** (exponent - 1), modulus)
    return tuple(pow(root, a, modulus) for a in range(1, h) if gcd(a, h) == 1)


class Norms:
    """Recover signed integer norms by evaluations in a large splitting field."""

    def __init__(self, n):
        orders = set(divisors(n) + divisors(n - 1))
        self.degree = {d: fmpz_poly.cyclotomic(d).degree() for d in orders}
        # Also handle each of the three halves in a fully paired partition.
        # One such half can have h-2 terms, more than the smallest full part.
        bound = 2 * (n // 2 + 1) ** max(self.degree.values())
        step = n * (n - 1)
        modulus = 1 + (bound // step + 1) * step
        while not fmpz(modulus).is_prime():
            modulus += step
        self.modulus = modulus
        self.half = n // 2
        self.low_mask = (1 << self.half) - 1
        self.tables = {}
        for d in orders:
            tables = []
            for root in roots_of_order(d, modulus, 1):
                powers = [pow(root, i, modulus) for i in range(n)]
                tables.append(
                    (
                        half_sums(powers[: self.half], modulus),
                        half_sums(powers[self.half :], modulus),
                    )
                )
            self.tables[d] = tables

    def at(self, d, mask):
        lo, hi = mask & self.low_mask, mask >> self.half
        f = g = 1
        for a, b in self.tables[d]:
            value = a[lo] + b[hi]
            f = f * value % self.modulus
            g = g * (value - 1) % self.modulus
        return min(f, self.modulus - f), min(g, self.modulus - g)


@lru_cache(None)
def compatible(d, value):
    """Keep only primes at which Phi_d can have an integer root."""
    result = 1
    for prime, exponent in factorization(value):
        s = int(prime)
        h = d
        while h % s == 0:
            h //= s
        if (s - 1) % h == 0:
            result *= s ** int(exponent)
    return result


@lru_cache(None)
def lower_ratio(orders):
    """Certify H(Q) >= Q^deg(H)/ratio for Q >= 2 by shifting the polynomial."""
    poly = fmpz_poly([1])
    for d in orders:
        poly *= fmpz_poly.cyclotomic(d)
    x = fmpz_poly([2, 1])
    ratio = 1
    while any(c < 0 for c in (ratio * poly(x) - x ** poly.degree()).coeffs()):
        ratio *= 2
    return ratio


@lru_cache(None)
def check_remaining_weight(n, q):
    if len(factorization(q)) == 1:
        raise ArithmeticError(("unresolved prime-power base", n, q))


SMALL_PRIMES = tuple(p for p in range(2, 501) if fmpz(p).is_prime())


def passes_carries(weight, parts, upper, assigned, cap):
    for r in SMALL_PRIMES:
        if upper % r and assigned % r:
            continue
        power = r
        while power <= weight:
            if weight // power != sum(v // power for v in parts):
                while upper % r == 0:
                    upper //= r
                while assigned % r == 0:
                    assigned //= r
                break
            power *= r
    lower = (weight - 1) // assigned
    return min(upper, cap) ** 3 >= 27 * lower**2


def evaluate(mask, powers):
    total = 0
    while mask:
        bit = mask & -mask
        total += powers[bit.bit_length() - 1]
        mask -= bit
    return total


def unpaired(n, norms):
    h = n // 2
    # Bound the norm at order 2h of any part with nonzero norm there.
    norm_cap = (h * h) ** ((h - 1) // 2) // (h - 1) ** ((h - 1) // 2)
    count = largest = tests = 0
    for k in range(2, n // 3 + 1):
        for positions in combinations(range(n), k):
            mask = sum(1 << i for i in positions)
            c = compatible(h, norms.at(h, mask)[0])
            norm, _ = norms.at(n, mask)
            c *= compatible(n, norm) if norm else norm_cap
            linear, _ = norms.at(2, mask)
            z = int(linear == 0)
            if linear:
                c *= linear
            count += 1
            if not z and c < 3:
                continue
            orders = []
            fixed = 1
            for d in divisors(n - 1):
                a, b = norms.at(d, mask)
                if a and b:
                    orders.append(d)
                    fixed *= compatible(d, a) * compatible(d, b)
            e = sum(norms.degree[d] for d in orders)
            assert 2 * e > 3 * z
            # (Q+1)/Q <= 3/2 for Q >= 2.
            bound = (c * 3**z) ** 3 * (lower_ratio(tuple(orders)) * fixed) ** 2 // (27 * 8**z)
            limit = floor_root(bound, 2 * e - 3 * z)
            largest = max(largest, limit)
            for q in range(2, limit + 1):
                weight, powers = powers_at(n, q)
                f = evaluate(mask, powers)
                u = min(gcd(f, weight), c * (q + 1) ** z)
                t = (weight - 1) // gcd(weight - 1, f * (f - 1))
                tests += 1
                if (
                    u >= 3
                    and 27 * t * t <= u**3
                    and passes_carries(
                        weight, (f, weight - f), gcd(f, weight), (weight - 1) // t, u
                    )
                ):
                    check_remaining_weight(n, q)
    assert count == sum(comb(n, k) for k in range(2, n // 3 + 1))
    return count, largest, tests


@lru_cache(maxsize=100000)
def zeros(mask, d, s, j):
    modulus = s**j
    return frozenset(
        r
        for r in roots_of_order(d, s, j)
        if sum(pow(r, i, modulus) for i in range(mask.bit_length()) if mask & (1 << i)) % modulus
        == 0
    )


def choices(d, cap, parts, required):
    for p, v in factorization(cap):
        s = int(p)
        options = [(1, 0)]
        for j in range(1, int(v) + 1):
            a, b, c = (zeros(mask, d, s, j) for mask in parts)
            roots = a & b & c if required == 3 else (a & b) | (a & c) | (b & c)
            if not roots:
                break
            options.extend((s**j, r) for r in roots)
        if len(options) > 1:
            yield options


def combine(classes):
    modulus, residue = 1, 0
    for w, r in classes:
        common = gcd(modulus, w)
        if (r - residue) % common:
            return None
        v = w // common
        shift = (r - residue) // common * pow(modulus // common, -1, v) % v
        residue += modulus * shift
        modulus *= v
    return modulus, residue


def paired_limit(h, u, c, degree_gap):
    # (Q^h+1)/Q^h <= (2^h+1)/2^h.
    bound = (2**h + 1) ** 3 * u**3 * c**2 // (27 * 2 ** (3 * h))
    return floor_root(bound, degree_gap)


def base_classes(h, parts, u, orders, constants, gap):
    scale = lower_ratio(orders)
    limit = paired_limit(h, u, scale * prod(constants), gap)
    if limit <= 200:
        yield 1, 0, limit
        return
    weights = list(choices(h, u, parts, 3))
    shifts = [options for d, c in zip(orders, constants) for options in choices(d, c, parts, 2)]
    for chosen in product(*weights, *shifts):
        congruence = combine(chosen)
        if congruence is None:
            continue
        u1 = prod(w for w, _ in chosen[: len(weights)])
        c1 = scale * prod(w for w, _ in chosen[len(weights) :])
        yield congruence[0], congruence[1], paired_limit(h, u1, c1, gap)


def paired(n, norms):
    h = n // 2
    full = (1 << h) - 1
    ds = divisors(n - 1)
    data = {}
    for mask in range(1, full):
        weight = compatible(h, norms.at(h, mask)[0])
        bounds = []
        for d in ds:
            value, _ = norms.at(d, mask)
            bounds.append(compatible(d, value) if value else 0)
        data[mask] = weight, bounds

    count = largest = tests = 0
    for x in range(1, full):
        rest = full ^ x
        y = rest
        while y:
            y = (y - 1) & rest
            z = rest ^ y
            if not x < y < z:
                continue
            count += 1
            parts = (x, y, z)
            records = [data[i] for i in parts]
            u = gcd(*(a for a, _ in records))
            orders = []
            constants = []
            for d, (a, b, c) in zip(ds, zip(*(b for _, b in records))):
                bound = lcm(gcd(a, b), gcd(b, c), gcd(c, a))
                if bound:
                    orders.append(d)
                    constants.append(bound)
            gap = 2 * sum(norms.degree[d] for d in orders) - 3 * h
            assert gap > 0
            for modulus, residue, limit in base_classes(
                h, parts, u, tuple(orders), constants, gap
            ):
                largest = max(largest, limit)
                start = 2 + (residue - 2) % modulus
                for q in range(start, limit + 1, modulus):
                    weight, powers = powers_at(n, q)
                    a, b, c = (evaluate(f, powers) for f in parts)
                    upper = (1 + q**h) * gcd(a, b, c)
                    lower = (weight - 1) // gcd(weight - 1, lcm(gcd(a, b), gcd(b, c), gcd(c, a)))
                    tests += 1
                    if 27 * lower**2 <= upper**3 and passes_carries(
                        weight,
                        tuple((1 + q**h) * v for v in (a, b, c)),
                        upper,
                        (weight - 1) // lower,
                        upper,
                    ):
                        check_remaining_weight(n, q)
    assert count == (3**h - 3 * 2**h + 3) // 6
    return count, largest, tests


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("lengths", nargs="*", type=int, default=[10, 14, 22, 26])
    parser.add_argument("--partitions", choices=("both", "paired", "unpaired"), default="both")
    args = parser.parse_args()
    for n in args.lengths:
        assert n % 2 == 0 and n // 2 >= 3 and fmpz(n // 2).is_prime() and n // 2 % 2
        norms = Norms(n)
        for routine in (unpaired, paired):
            if args.partitions in ("both", routine.__name__):
                print(n, routine.__name__, *routine(n, norms), flush=True)
                check_remaining_weight.cache_clear()
