"""Certify fixed column ratios with independent column and layer spacings."""
from argparse import ArgumentParser
from collections import Counter
from math import gcd, prod

from check_geometric_product import MULTIPLICITIES, roots, sectors
from check_seven_span import companion, coprime_part, odd_part
from check_translated_triple_columns import allocations
from check_translated_triple_lines import specialize

PRIMES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61)


def check(s, t):
    assert 0 < s < t and gcd(s, t) == 1
    grid = tuple((c, layer) for layer in range(2) for c in (0, s, t))
    counts = Counter()
    observed_roots = set()

    def height(v):
        return 2*(t*v[0]+v[1])

    for xs, ys, zi in allocations():
        x, y = (tuple(grid[i] for i in indices) for indices in (xs, ys))
        z = grid[zi]
        pairs = [(x, y), (x, tuple(sorted(y+(z,)))), (y, tuple(sorted(x+(z,))))]
        pending = set()
        for v, w, fs, degree in sectors((1, 0), (0, 1), pairs, height):
            counts['oriented sectors'] += 1
            positive = {root for f in fs for root in roots(f)}
            observed_roots.update(positive)
            if 1 in positive and all(v):
                for p in PRIMES:
                    b = 1
                    while p**b < 64:
                        pending.add((p, b*v[0], b*v[1]))
                        b += 1
            for p in PRIMES:
                for value in positive - {1}:
                    ell, power = 0, 1
                    while power < value:
                        ell += 1
                        power *= p
                    if power != value:
                        continue
                    r = ell*(t*w[0]+w[1])
                    e = (ell*w[0], ell*w[1])
                    f, g, m_poly = (specialize(p, block, v, e) for block in (x,y,grid))
                    h = m_poly.degree()
                    margin = 2*h - 3*f.gcd(g).degree() - 2*m_poly.gcd(f*g).degree()
                    assert margin > 0
                    limit = h*(30+16*r)//margin
                    counts['nonhomogeneous root checks'] += 1
                    for n in range(ell, limit+1):
                        pending.add((p, n*v[0]+ell*w[0], n*v[1]+ell*w[1]))
            norm = prod(sum(map(abs, f))**j for f, j in zip(fs, MULTIPLICITIES))
            gap = height(v)+height(w)-degree
            assert gap > 0
            budget = (norm//27).bit_length()
            for b in range(budget//height(v)+1):
                for ell in range((budget-b*height(v))//gap+1):
                    d, k = (b+ell)*v[0]+ell*w[0], (b+ell)*v[1]+ell*w[1]
                    if not d or not k:
                        continue
                    exponent = b*height(v)+ell*gap
                    for p in PRIMES:
                        if 27*p**exponent > norm:
                            break
                        pending.add((p, d, k))
        for p, d, k in pending:
            if (p == 2 and k in (s*d,t*d,(t-s)*d)) or p**gcd(d,k) >= 64:
                continue
            a, b, c = p**(s*d), p**(t*d), p**k
            values = (1,a,b,c,a*c,b*c)
            m = sum(values)
            x0, y0 = (sum(values[i] for i in indices) for indices in (xs, ys))
            u = coprime_part(odd_part(gcd(x0,y0)), p*m)
            cap = odd_part(gcd(m, x0*y0))
            counts['integer checks'] += 1
            if u**3 * cap**2 >= 27*m*m:
                counts['companion tests'] += 1
                companion(m, u, cap, counts)
    print((s,t), dict(counts), 'positive roots', sorted(observed_roots), flush=True)


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--max-column', type=int, default=5)
    args = parser.parse_args()
    for t in range(3, args.max_column+1):
        for s in range(1, t):
            if gcd(s,t) == 1:
                check(s,t)
