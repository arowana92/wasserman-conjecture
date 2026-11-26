"""Test rational-line divisibilities and primitive displacement bounds with exact coefficients."""

from collections import Counter
from itertools import combinations, combinations_with_replacement, permutations, product
from math import comb, gcd, lcm, prod


def cross(x, y):
    return tuple(x[(i+1)%3]*y[(i+2)%3]-x[(i+2)%3]*y[(i+1)%3] for i in range(3))


def primitive(values):
    factor = gcd(*values)
    return tuple(v//factor for v in values)


def check_line(n, normal, coefficient):
    den = gcd(normal[1]-normal[0], normal[2]-normal[0])
    assert den > 0 and n % den == 0 and coefficient % den == 0
    if 0 not in normal:
        L = lcm(*map(abs, normal))
        assert coefficient % (den*((n-1)//gcd(n-1, L))) == 0


def ratios(row):
    result = []
    for a, b in combinations(row, 2):
        g = gcd(a, b)
        result.append((a//g, b//g, g))
    return result


def verify(limit=40):
    counts = Counter()
    min_shared = None
    for n in range(3, limit+1):
        rows = [(a,b,n-a-b) for a in range(1,n//3+1)
                for b in range(a,(n-a)//2+1)]
        coefficients = {r: comb(n,r[0])*comb(n-r[0],r[1]) for r in rows}
        for A in rows:
            g = gcd(n*(n-1), coefficients[A])
            assert g**3 >= 27*(n-1)**2
            # Independent line tests, including many signed and zero coordinates.
            for x, y in product(range(-3,4), repeat=2):
                normal = cross(A,(x,y,1))
                if any(normal):
                    check_line(n, primitive(normal), coefficients[A])
        for A,C in combinations_with_replacement(rows,2):
            KA, KC = coefficients[A], coefficients[C]
            D, J = gcd(n, KA, KC), gcd(n-1, KA, KC)
            G = D*J
            for (a,b,u),(c,d,v) in product(ratios(A),ratios(C)):
                L = lcm(a,b,c,d)
                assert (G*u*v*L) % (n*(n-1)) == 0
                if a == c and b == d:
                    assert G > 4
                    min_shared = G if min_shared is None else min(G,min_shared)
                    counts['shared-ratio comparisons'] += 1
                counts['ratio comparisons'] += 1
            for other in set(permutations(C)):
                if A == other:
                    continue
                raw = cross(A, other)
                normal = primitive(raw)
                assert sum(a*v for a,v in zip(A,normal)) == 0
                assert sum(c*v for c,v in zip(other,normal)) == 0
                check_line(n, normal, G)
                delta = tuple(a-b for a,b in zip(A,other))
                height = max(map(abs,primitive(delta)))
                den = gcd(normal[1]-normal[0], normal[2]-normal[0])
                assert max(normal)-min(normal) == den*height
                if 0 not in normal:
                    assert 27*(n-1) <= 4*J*D**3*height**3
                    counts['nonzero normal comparisons'] += 1
                else:
                    assert G > 4
                counts['aligned row comparisons'] += 1
    print(dict(counts), 'least endpoint gcd for a shared ratio:', min_shared)


if __name__ == '__main__':
    verify()
