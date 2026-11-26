"""Check the cyclotomic factor and the extra offsets in bases two and three."""
from functools import cache
from math import gcd, lcm
from random import Random

from flint import fmpz_poly


def doubled_offsets(p):
    gaps = (1, 2) if p == 2 else (1,) if p == 3 else ()
    return [(0, 0, 0)] + [v for g in gaps for v in ((0,g,0), (0,0,g), (0,-g,-g))]


def strip(terms):
    low = min(i for i, c in terms)
    coeffs = [0] * (max(i for i, c in terms)-low+1)
    for i, c in terms:
        coeffs[i-low] += c
    return fmpz_poly(coeffs)


@cache
def line_factor(p, v, w):
    def build(terms):
        low = min(a for _, a in terms)
        return strip([(i, p**(a-low)) for i, a in terms])
    f = build(((0,0),(v[0],w[0])))
    g = build(((0,0),(v[1],w[1]),(v[2],w[2])))
    common = f.gcd(g)
    common = common // gcd(*(int(c) for c in common.coeffs()))
    assert common.degree() <= 2
    return common


def principal(p, d, e, f):
    values = []
    for offset in doubled_offsets(p):
        delta = tuple(2*a-h for a, h in zip((d,e,f), offset))
        g = gcd(*delta)
        v = tuple(a//g for a in delta)
        b = g//2
        w = tuple(a-b*i for a, i in zip((d,e,f), v))
        poly = line_factor(p,v,w)
        value = abs(int(poly(p**b)))
        while value % p == 0:
            value //= p
        values.append(value)
    return lcm(*values), values[0]


def cyclotomic(p, d, e, f):
    common = gcd(d,e,f)
    if (d//common) % 6 == 3 and {(e//common) % 6,(f//common) % 6} == {2,4}:
        t = p**common
        return t*t-t+1
    return 1


def check():
    cases = extras = 0
    witnesses = {}
    for p in (2,3,5,7,11):
        powers = [p**i for i in range(25)]
        for d in range(1,25):
            for f in range(25):
                for e in range(f+1):
                    v, homogeneous = principal(p,d,e,f)
                    g = gcd(1+powers[d],1+powers[e]+powers[f])
                    assert homogeneous == cyclotomic(p,d,e,f)
                    assert g % v == 0
                    if p >= 5:
                        assert v == homogeneous
                    if v > homogeneous:
                        extras += 1
                        if 0 < e < f and v >= 100:
                            witnesses.setdefault(p,(d,e,f,v,g))
                    cases += 1
    print(f"{cases} exact comparisons, {extras} strict increases from small-prime offsets")
    print("Distinct positive examples", witnesses)
    rng = Random(785)
    for case in range(400):
        p = (5,7,11,13)[case % 4]
        if case % 3:
            d = rng.randint(8000, 12000)
            e, f = sorted(rng.sample(range(1,12001), 2))
        else:
            b = rng.randint(1800, 3000)
            d, e, f = 3*b, 2*b, 4*b
        h = max(d,f)
        scale = 1
        while (scale+1)**3 <= h*h:
            scale += 1
        v = cyclotomic(p,d,e,f)
        g = gcd(1+p**d, 1+p**e+p**f)
        assert g % v == 0 and g//v <= p**(13*scale)
    for b in (11,91,501,5001):
        for p, d, e, f in ((2,b,b+1,2*b), (3,3*b,(b+1)//2,b)):
            g = gcd(1+p**d, 1+p**e+p**f)
            assert g == (1+2**b if p == 2 else 1+p**e+p**f)
    print("400 large-span comparisons and eight exact small-prime family checks")


if __name__ == '__main__':
    check()
