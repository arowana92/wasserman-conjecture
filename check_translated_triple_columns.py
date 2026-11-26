"""Certify a bounded normalized column span at arbitrary shifts and translations."""
from argparse import ArgumentParser
from collections import Counter
from itertools import combinations
from math import gcd, isqrt

from check_seven_span import companion, coprime_part, odd_part


def allocations():
    """Generate the 36 surviving labelled choices from the three row types."""
    for z in range(6):
        for x in combinations((s for s in range(6) if s != z), 2):
            if x[0] % 3 == x[1] % 3:
                continue
            # Reject X+Y and Z(X+Y) in the coordinates relative to z.
            if x[0] // 3 == x[1] // 3 and all(s % 3 != z % 3 for s in x):
                continue
            y = tuple(s for s in range(6) if s != z and s not in x)
            yield x, y, z


def coefficients(terms, columns):
    return tuple(sum(columns[s % 3] for s in terms if s // 3 == layer)
                 for layer in range(2))


def linear_cap(f, g):
    """A fixed multiple of the gcd away from the evaluation base, or zero."""
    constants = [a+b for a, b in (f, g) if a*b == 0]
    if constants:
        return gcd(*constants)
    return abs(f[0]*g[1]-f[1]*g[0])


def check(base, p, width):
    counts = Counter()
    powers = [base**j for j in range(6*width+9)]
    layouts = tuple(allocations())
    assert len(layouts) == 36
    for j in range(2, width+1):
        lo = max(1, 101-j)
        hi = 6*j+7
        # Below the span-100 theorem or beyond the seven-unit gap bound.
        if lo > hi:
            continue
        for i in range(1, j):
            columns = (1, powers[i], powers[j])
            total = sum(columns)
            for xx, yy, z_index in layouts:
                x = coefficients(xx, columns)
                y = coefficients(yy, columns)
                z = coefficients((z_index,), columns)
                caps = [linear_cap(x, y),
                        linear_cap(x, tuple(a+b for a, b in zip(y, z))),
                        linear_cap(y, tuple(a+b for a, b in zip(x, z)))]
                # Every determinant is bounded by 9*base**(2*j). Strip its
                # base prime and 2 since none can divide an assigned factor.
                caps = [odd_part(c // gcd(c, powers[2*j+8])) if c else 0 for c in caps]
                if caps[0] == 1:
                    continue
                if all(caps):
                    u, a, c = caps
                    maximum = isqrt(u**3*(a*c)**2//27)//total
                else:
                    maximum = None
                for k in range(lo, hi+1):
                    t = powers[k]
                    if maximum is not None and 1+t > maximum:
                        break
                    if p == 2 and k in (i, j, j-i):
                        continue
                    bx, by = x[0]+x[1]*t, y[0]+y[1]*t
                    m = total*(1+t)
                    g = coprime_part(odd_part(gcd(caps[0], bx, by)), m*p)
                    if g < (3 if p == 2 else 7):
                        continue
                    s = odd_part(gcd(m, bx*by))
                    counts['numerical size tests'] += 1
                    if g**3*s*s >= 27*m*m:
                        counts['companion tests'] += 1
                        companion(m, g, s, counts)
    return counts


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--width', type=int, default=48)
    args = parser.parse_args()
    if args.width < 2:
        parser.error("the column span must be at least two")
    totals = Counter()
    for p in range(2, 64):
        if any(p % q == 0 for q in range(2, isqrt(p)+1)):
            continue
        base = p
        while base < 64:
            totals.update(check(base, p, args.width))
            totals['bases'] += 1
            base *= p
    print(dict(totals))
