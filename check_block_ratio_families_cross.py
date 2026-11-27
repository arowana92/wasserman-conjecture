"""Independent check by odd cofactors, pair sums, and possible single coordinates."""
from itertools import combinations
from math import gcd, prod


def powers(n):
    result = []
    candidate = 2
    while candidate*candidate <= n:
        if n % candidate == 0:
            w = 1
            while n % candidate == 0:
                n //= candidate
                w *= candidate
            result.append((candidate, w))
        candidate = 3 if candidate == 2 else candidate+2
    if n > 1:
        result.append((n, n))
    return result


def companion_exists(cap, atoms):
    t = prod(atoms)
    full = (1 << len(atoms))-1
    allowed = {}
    for mask in range(full+1):
        w = prod(q for j, q in enumerate(atoms) if mask & (1 << j))
        zero = t//w
        # The coordinate is zero modulo 'zero' and agrees with cap modulo w.
        least = zero*((cap*pow(zero, -1, w)) % w) if w > 1 else 0
        least = least or t
        if least+2 <= cap:
            allowed[mask] = least
    for a, x in allowed.items():
        for b, y in allowed.items():
            if a & b:
                continue
            c = full ^ (a | b)
            if c in allowed and x+y+allowed[c] <= cap:
                return True
    return False


def primitive_totals():
    pairs = set()
    for width in range(2, 7):
        for size in range(2, 7):
            positions = range(0, width*size, width)
            maximal = sum(2**j for j in positions)
            pairs.add((maximal, 2**width-1))
            if width+size <= 10:
                for count in range(size):
                    for extras in combinations(list(positions)[1:], count):
                        upper = 1+sum(2**j for j in extras)
                        for low in range(1, 2**width):
                            pairs.add((upper, low))
    values = set()
    for b, c in pairs:
        # Once both strict inequalities hold, all larger shifts are covered.
        h = 1
        while True:
            t = (b << h)+c
            if (1 << h)>c and 27*t*t > b**7*c**7:
                break
            values.add(t//(t & -t))
            h += 1
    return len(pairs), values


def verify(total, stats):
    square_bound = 27*total*total
    factorization = None
    union = (total-1) & total
    while union:
        odd = union//(union & -union)
        g = gcd(odd, total)
        while g > 1:
            odd //= g
            g = gcd(odd, total)
        if odd >= 3 and union & (union-1):
            # Fix the lowest occupied bit in y, counting each unordered pair once.
            free = union ^ (union & -union)
            x = free
            while x:
                y = union-x
                cap = gcd(x, odd)
                if cap >= 3:
                    stats[0] += 1
                    q = gcd(total, x)
                    assigned = q*gcd(total//q, y)
                    if cap**3*assigned**2 >= square_bound:
                        stats[1] += 1
                        if factorization is None:
                            factorization = powers(total)
                        forced = [w for p,w in factorization if
                                  (x % w and y % w) or
                                  (x % w and x % p == 0) or
                                  (y % w and y % p == 0)]
                        if 27*prod(forced)**2 <= cap**3 and all(2*w < cap for w in forced):
                            stats[2] += 1
                            assert not companion_exists(cap, forced), (total,x,y,cap,forced)
                x = (x-1) & free
        union = (union-1) & total


if __name__ == "__main__":
    families, values = primitive_totals()
    stats = [0, 0, 0]
    for n in values:
        verify(n, stats)
    print(families, "families,", len(values), "distinct odd cofactors,", stats)
    assert families == 1674 and len(values) == 39304
    assert stats == [8933277, 910, 558]
