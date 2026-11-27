"""Exhaust the bounded-height branch for carry-free three-power translates."""

from argparse import ArgumentParser
from collections import Counter, defaultdict
from math import gcd, isqrt, prod


VERIFIED_THROUGH = 10**9


def allocations(n, p):
    """Ordered carry-free distributions of the digits of n among three parts."""
    if n == 0:
        yield 0, 0, 0
        return
    upper, digit = divmod(n, p)
    for x, y, z in allocations(upper, p):
        for i in range(digit + 1):
            for j in range(digit - i + 1):
                yield p*x + i, p*y + j, p*z + digit - i - j


def p_free(x, p):
    assert x > 0
    if p == 2:
        return x // (x & -x)
    while x % p == 0:
        x //= p
    return x


def height(offset, p):
    power, digits = 1, 0
    while power < 2*offset + 2:
        power *= p
        digits += 1
    return 90*digits + 9


def unitary(n, divisor):
    """Product of the full prime powers of n dividing divisor."""
    candidate = gcd(n, divisor)
    missing = n // candidate
    while (shared := gcd(candidate, missing)) > 1:
        candidate //= shared
    return candidate


def inspect(n, row, caps, p):
    product = prod(caps[1:])
    target = 27*(n-1)**2
    if caps[0]**3 * product**2 < target:
        return None
    parts = [p_free(x, p) for x in row]
    if len(set(parts)) == 1:
        return "unshifted shapes"
    for j in range(3):
        if parts[(j+1)%3] == parts[(j+2)%3] == p_free(row[j]-1, p):
            return "shifted shapes"
    u, s = unitary(n, caps[0]), unitary(n-1, product)
    if u >= 3 and n//u >= 3 and u**3*s**2 >= target:
        return "unresolved rows"
    return None


def scan(p, offsets):
    groups = defaultdict(list)
    for offset in offsets:
        if offset*(offset-1) % p == 0:
            for x, y, z in allocations(offset, p):
                groups[y, z].append((offset, x, height(offset, p)))
    counts = Counter()
    unresolved = []
    if not groups:
        return counts, unresolved

    limit = max(offsets)
    assert (2*limit)**5 < 27*VERIFIED_THROUGH**2
    begin = 0
    while (p**begin + limit)**5 <= 27*VERIFIED_THROUGH**2:
        begin += 1
    max_height = max(h for params in groups.values() for _, _, h in params)
    powers = [p**i for i in range(max_height + 1)]
    assert powers[begin] > limit
    strict = int(p == 2)
    for (y, z), params in groups.items():
        hmax = max(h for _, _, h in params)
        pmax = max(t // p_free(t, p) for c, _, _ in params for t in (c*(c-1),))
        # Precompute the prime-to-p parts of the two smaller entries and their shifts.
        small = [[p_free(power + delta, p) for power in powers[begin:hmax+1]]
                 for delta in (y, y-1, z, z-1)]
        for b in range(begin+strict, hmax+1-strict):
            f, fm = small[0][b-begin], small[1][b-begin]
            lower = max(VERIFIED_THROUGH, powers[b+strict])
            target_bits = (27*lower**2).bit_length()
            for c in range(begin, b+1-strict):
                t, tm = small[2][c-begin], small[3][c-begin]
                g, left, right = gcd(f, t), gcd(fm, t), gcd(f, tm)
                if (pmax**3).bit_length() + 3*g.bit_length() + 2*(left.bit_length()+right.bit_length()) < target_bits:
                    continue
                upper = 1 + isqrt(pmax**3 * g**3 * (left*right)**2 // 27)
                subtotal = powers[b] + powers[c]
                if upper < max(VERIFIED_THROUGH+1, subtotal+powers[b+strict]+2):
                    continue
                middle, last = powers[b]+y, powers[c]+z
                both = gcd(middle, last)
                minus_middle, minus_last = gcd(middle-1, last), gcd(middle, last-1)
                for offset, x, top in params:
                    counts["pair candidates"] += 1
                    for a in range(b+strict, top+1):
                        n = subtotal + powers[a] + offset
                        if n > upper:
                            break
                        if n <= VERIFIED_THROUGH or (p == 3 and a == b == c):
                            continue
                        counts["rows"] += 1
                        first = powers[a] + x
                        caps = (gcd(first,both), gcd(first-1,both),
                                gcd(first,minus_middle), gcd(first,minus_last))
                        result = inspect(n, (first,middle,last), caps, p)
                        if result:
                            counts[result] += 1
                        if result == "unresolved rows":
                            unresolved.append((p, offset, a, b, c, x, y, z))
    return counts, unresolved


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--limit", type=int, default=63)
    parser.add_argument("--one", type=int, help="check a single offset")
    args = parser.parse_args()
    offsets = [args.one] if args.one is not None else list(range(2,args.limit+1))
    if not offsets or min(offsets) < 2 or max(offsets) > 2048:
        parser.error("choose offsets in 2..2048")
    totals = Counter()
    for p in range(2,max(offsets)+1):
        if any(p % d == 0 for d in range(2,isqrt(p)+1)):
            continue
        counts, pending = scan(p, offsets)
        if counts:
            print(p, dict(counts), flush=True)
        for example in pending:
            print("RETAINED", example, flush=True)
        totals.update(counts)
    print("total", dict(totals))
