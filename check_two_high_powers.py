"""Exhaust the finite range for offsets 3 through 8 and arbitrary outer exponent."""

from math import gcd

from prime_certificates import primes_through


def allocations(c):
    # The two tags index (0, H, J, H+J). The other entries are offsets on x,y.
    for a in range(c+1):
        for b in range(c-a+1):
            z = c-a-b
            if min(a, b, z) > 0:
                yield 1, 2, a, b
            # If one unshifted slot is constant put it first. The common
            # scale divides this slot, so its value must be at least three.
            if a < 3:
                continue
            if a <= b:
                yield 0, 0, a, b
            if b > 0:
                yield 0, 1, a, b
                yield 0, 2, a, b
            if z > 0:
                yield 0, 3, a, b


def check(c):
    largest_r = (2*c)**5
    largest_h = (2*c)**64
    choices = tuple(allocations(c))
    weights = 0
    survivors = []
    for p in primes_through(largest_r):
        if p <= c:
            continue
        r = p
        while r <= largest_r:
            j, high = 2, r*r
            while high <= largest_h:
                low = 1
                for i in range(1, j):
                    low *= r
                    if gcd(i, j) != 1:
                        continue
                    weights += 1
                    d = high+low+c
                    units = 0, high, low, high+low
                    for left, middle, a, b in choices:
                        x = units[left]+a
                        y = units[middle]+b
                        common = gcd(x, y)
                        if common < 3:
                            continue
                        assigned = gcd(d, x*y)
                        if common**3 * assigned**2 < 27*d*d:
                            continue
                        remove = gcd(common, p*d)
                        while remove > 1:
                            common //= remove
                            remove = gcd(common, p*d)
                        if common**3 * assigned**2 >= 27*d*d:
                            survivors.append((p, r, j, i, d, x, y, common, assigned))
                j += 1
                high *= r
            r *= p
    return weights, survivors


def low_digits_fit(x, y, z, p):
    while x or y or z:
        if x % p + y % p + z % p >= p:
            return False
        x //= p
        y //= p
        z //= p
    return True


def small_base_check(c):
    checked = 0
    for p in primes_through(c):
        choices = [entry for entry in allocations(c)
                   if low_digits_fit(entry[2], entry[3], c-entry[2]-entry[3], p)]
        high = p*p
        while high <= (2*c)**64:
            low = p
            while high > c*low:
                d = high+low+c
                if low <= c:
                    # The single-high-power theorem already excludes larger d.
                    if 27*d*d >= (low+c)**5:
                        low *= p
                        continue
                    pairs = ((x, y) for x in range(1, d-1) for y in range(x, d-x)
                             if low_digits_fit(x, y, d-x-y, p))
                else:
                    units = (0, high, low, high+low)
                    pairs = ((units[i]+a, units[j]+b) for i, j, a, b in choices)
                for x, y in pairs:
                    scale = gcd(x, y)
                    assigned = gcd(d, x*y)
                    if scale**3 * assigned**2 >= 27*d*d:
                        while gcd(scale, p*d) > 1:
                            scale //= gcd(scale, p*d)
                        assert scale**3 * assigned**2 < 27*d*d, (p, c, high, low, x, y)
                checked += 1
                low *= p
            high *= p
    return checked


if __name__ == '__main__':
    for offset in range(3, 9):
        weights, survivors = check(offset)
        expected = [(13, 13, 2, 1, 190, 171, 18, 9, 38)] if offset == 8 else []
        assert survivors == expected, (offset, survivors)
        for p, _, _, _, d, _, _, common, _ in survivors:
            assert common == 9 and p % 3 == d % 3 == 1
        small = small_base_check(offset)
        print(f'Offset {offset}: {weights} large-base and {small} small-base inner pairs, '
              f'{len(survivors)} size-test survivors, zero surviving assignments.', flush=True)
