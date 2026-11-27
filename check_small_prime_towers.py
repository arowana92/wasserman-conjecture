"""Finite periods for the remaining small-prime exponent families."""

from functools import lru_cache
from math import gcd, lcm

# d, prime bases, first slopes, companion slopes, W residue and modulus,
# and the prime bases that the first row must avoid (shifted_directions.md).
CASES = (
    (380, (1733, 1973, 5717), (130, 247, 3), (1, 4, 8), 823, 2470, (5, 13, 19)),
    (629, (4943,), (259, 357, 13), (1, 2, 4), 5419, 8806, (7, 17, 37)),
    (780, (6379,), (224, 455, 101), (1, 3, 3), 919, 1820, (2, 5, 7, 13)),
    (850, (2851, 3271, 8731, 8971), (175, 476, 199), (1, 2, 4), 5651, 5950, (5, 7, 17)),
    (880, (3889,), (272, 561, 47), (2, 5, 10), 2801, 2992, (2, 11, 17)),
    (950, (4463,), (275, 627, 48), (1, 2, 8), 4463, 10450, (5, 11, 19)),
    (962, (8861,), (143, 814, 5), (1, 2, 8), 6349, 10582, (11, 13, 37)),
    (986, (2447,), (119, 203, 664), (1, 2, 4), 4215, 6902, (7, 17, 29)),
    (986, (2447,), (203, 714, 69), (1, 2, 4), 3501, 6902, (7, 17, 29)),
)


@lru_cache(maxsize=256)
def orbit(q, modulus):
    assert gcd(q, modulus) == 1
    value = 1 % modulus
    powers = [value]
    while True:
        value = value * q % modulus
        if value == powers[0]:
            return powers
        powers.append(value)


def free(parts, p, modulus):
    while modulus > 1:
        if sum(x % modulus for x in parts) >= modulus:
            return False
        modulus //= p
    return True


def local(case, p, w, W):
    d, _, (a, b, c), slopes, target, _, forced = case
    if p == 0:
        return W == target
    n = 1 + d * W
    U = sum(slopes)
    factor = gcd(U, w)
    reduced = w // factor
    assert reduced % p == 0 and gcd(U // factor, reduced) == 1
    if n % factor:
        return False
    V = (n // factor) * pow(U // factor, -1, reduced) % reduced
    if U % p == 0 and V % p == 0:
        return False
    first = free((a * W, b * W, c * W + 1), p, w)
    other = free(tuple(f * V for f in slopes), p, reduced)
    first_required = p in forced or U % p == 0
    return ((not first_required or first) and (V % p != 0 or other)
            and (d % p != 0 or first or other))


def floor_check(case, p, w, W):
    # Independently choose an integer lift of W that makes the entire companion
    # integral. Compare factorial floors at the digits determined by the modulus.
    d, _, slopes, other, target, _, forced = case
    if p == 0:
        return W == target
    U = sum(other)
    for lift in range(W, W + w * U, w):
        n = 1 + d * lift
        if n % U == 0:
            break
    else:
        return False
    V = n // U
    if gcd(U, V) % p == 0:
        return False
    A = [s * lift for s in slopes]
    A[-1] += 1
    C = [f * V for f in other]
    bound = w // gcd(U, w)
    first_val = other_val = 0
    layer = p
    while layer <= w:
        first_val += n // layer - sum(x // layer for x in A)
        if layer <= bound:
            other_val += n // layer - sum(x // layer for x in C)
        layer *= p
    if first_val and (p in forced or U % p == 0):
        return False
    if other_val and V % p == 0:
        return False
    return d % p != 0 or min(first_val, other_val) == 0


def tests(case):
    yield 0, case[5]
    seen = set()
    for cap in (256, 4096, 65536):
        for p in (2, 3, 5, 7, 11, 13, 17, 19, 29, 37):
            w = p
            while w * p <= cap:
                w *= p
            if w // gcd(w, sum(case[3])) == 1:
                w *= p
            if w not in seen:
                seen.add(w)
                yield p, w


def allowed_indices(case, q, depth, p, w, multipliers):
    values = orbit(q, w)
    accepted = []
    for x in values:
        result = local(case, p, w, x)
        assert result == floor_check(case, p, w, x)
        accepted.append(result)
    stride = len(values)
    inner = orbit(q, stride)
    if depth == 2:
        inner = [inner[i] for i in orbit(q, len(inner))]
    # Check the periods' entries by independent modular exponentiation.
    assert all(x == pow(q, j, w) for j, x in enumerate(values))
    period = len(inner)
    for i, e in enumerate(inner):
        exponent = i if depth == 1 else pow(q, i, len(orbit(q, stride)))
        assert e == pow(q, exponent, stride)
    for b in multipliers:
        residues = {i for i, e in enumerate(inner) if accepted[b * e % stride]}
        yield b, (period, residues)


def clashes(first, second):
    m, a = first
    n, b = second
    g = gcd(m, n)
    return {i % g for i in a}.isdisjoint(i % g for i in b)


def sieve(case, q, depth, multipliers):
    choices = {b: {} for b in multipliers}
    for p, w in tests(case):
        if not choices:
            break
        for b, current in allowed_indices(case, q, depth, p, w, list(choices)):
            prior = choices[b]
            if not current[1] or any(clashes(current, old) for old in prior.values()):
                del choices[b]
            else:
                prior[w] = current
    return choices


def run():
    survivors = set()
    row_pairs = split = 0
    for case in CASES:
        d, bases, slopes, _, _, _, _ = case
        assert sum(slopes) == d
        for q in bases:
            row_pairs += 1
            # The existing cube theorem leaves only 986 when 3 divides the
            # exponent. The square theorem already treats every even multiplier.
            candidates = [b for b in range(1, 1001, 2) if d == 986 or b % 3]
            alive = sieve(case, q, 1, candidates)
            if d == 986 and 891 in alive:
                groups = [alive[891][w] for w in (4096, 2401, 24389)]
                period = lcm(*(m for m, _ in groups))
                assert period == 68208
                assert not any(all(i % m in allowed for m, allowed in groups)
                               for i in range(period))
                del alive[891]
            survivors.update((d, q, b) for b in alive)
            assert not sieve(case, q, 2, alive)
            split += len(alive)
    expected = {(629, 4943, 97), (629, 4943, 685), (986, 2447, 333), (986, 2447, 687)}
    assert survivors == expected and split == 4
    print(f"{row_pairs} row/base cases. Possible (d, q, b) at one level:", sorted(survivors))
    print("Every case excluded at two levels.")


if __name__ == "__main__":
    run()
