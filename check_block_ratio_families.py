"""Certify binary two-block families at every inner and outer exponent."""
from itertools import product
from math import gcd, prod

from check_endpoint_pair_bound import prime_powers


def submasks(value):
    mask = value
    while mask:
        yield mask
        mask = (mask-1) & value


def possible_companion(scale, powers):
    modulus = prod(powers)
    residues = [scale*(modulus//w)*pow(modulus//w, -1, w) % modulus for w in powers]
    for slots in product(range(3), repeat=len(powers)):
        parts = [0, 0, 0]
        for value, slot in zip(residues, slots):
            parts[slot] += value
        if sum((a-1) % modulus+1 for a in parts) <= scale:
            return True
    return False


def remaining_rows(total, counts):
    target = 27*total*total
    factors = None
    for x in submasks(total):
        if 2*x >= total:
            continue
        for y in submasks(total ^ x):
            if y < x or x+y == total:
                continue
            counts[0] += 1
            cap = gcd(x, y)
            cap //= cap & -cap
            if cap < 3:
                continue
            assigned = gcd(total, x*y)
            if cap**3*assigned**2 < target:
                continue
            common = gcd(cap, total)
            while common > 1:
                cap //= common
                common = gcd(cap, total)
            if cap**3*assigned**2 < target:
                continue
            counts[1] += 1
            if factors is None:
                factors = list(prime_powers(total))
            forced = [w for p, w in factors if not (
                (x % w == 0 and (y % w == 0 or y % p != 0)) or
                (y % w == 0 and x % p != 0))]
            if 27*prod(forced)**2 > cap**3 or any(2*w+1 > cap for w in forced):
                continue
            counts[2] += 1
            if possible_companion(cap, forced):
                yield x, y, cap, tuple(forced)


def separated(B, c):
    def distances(n):
        support = [i for i in range(n.bit_length()) if n >> i & 1]
        return {i-j for i in support for j in support if j<i}
    return distances(B).isdisjoint(distances(c))


def certify(B, c, cache, counts):
    assert separated(B,c)
    limit = (B*c)**7
    power = 2
    inner = 1
    residual = []
    while power <= c or 27*(B*power+c)**2 <= limit:
        total = B*power+c
        if total not in cache:
            cache[total] = tuple(remaining_rows(total, counts))
        if cache[total]:
            residual.append((inner,total,cache[total]))
        inner += 1
        power *= 2
    return residual


def families():
    small = set()
    for width in range(2, 7):
        for size in range(2, 7):
            if width+size > 10:
                continue
            upper = sum(1 << (i*width) for i in range(size))
            small.update((B,c) for B in submasks(upper) if B % 2
                         for c in range(1, 2**width))
    large = {(sum(1 << (i*width) for i in range(size)), 2**width-1)
             for width in range(2, 7) for size in range(2, 7)}
    return small, large


if __name__ == "__main__":
    small, large = families()
    cache = {}
    counts = [0, 0, 0]
    residual = []
    for B,c in sorted(small | large):
        failures = certify(B, c, cache, counts)
        if failures:
            residual.append((B,c,failures))
    print("Mask families",len(small | large),"distinct cofactors",len(cache),
          "ordered-by-size unshifted pairs",counts[0],flush=True)
    print("Passing product tests",counts[1],"passing individual power tests",counts[2])
    print("Unresolved families",len(residual))
    assert not residual, residual[:2]
    assert (len(small | large), len(cache), counts) == (1674, 71129, [92658606, 1351, 835])
