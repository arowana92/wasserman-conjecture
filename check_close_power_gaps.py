"""Certify the small offsets left by the two-high-power theorem."""

from collections import Counter
from itertools import product
from math import gcd

from check_six_digit_gaps import complement_by_parts, factors
from check_two_high_powers import low_digits_fit

CASES = ((5,3,1), (7,2,1), (7,2,2), (7,3,1), (7,5,1), (8,3,1), (8,5,1))
# Each entry is (b,v), for J+s = p^b (p^(h-b)+1) and t = p^v.
RAYS = {(5,3,1): {(0,0)}, (7,2,1): {(0,2),(1,0)}, (7,2,2): {(0,1)},
        (7,3,1): {(0,1)}, (7,5,1): {(0,0)}}
BIT_BOUND = 256


def determinant(x, y):
    return x[0]*y[1]-x[1]*y[0]


def split_forms(c, p, k):
    for high, low in product(range(3), repeat=2):
        slopes = tuple(k*(i == high)+(i == low) for i in range(3))
        for r in range(c+1):
            for s in range(c-r+1):
                t = c-r-s
                if not low_digits_fit(r, s, t, p):
                    continue
                x, y, z = tuple(zip(slopes, (r,s,t)))
                if x < y or (0,0) in (x,y,z):
                    continue
                # A constant unshifted slot is a multiple of U >= 3.
                # An unshifted pure p-power has no factor in common with U.
                if any((j == 0 and v < 3) or (j in (1,k) and v == 0)
                       for j, v in (x,y)):
                    continue
                yield x, y, z


def coprime_part(n, m):
    while gcd(n,m) > 1:
        n //= gcd(n,m)
    return n


def exclude_numeric(p, d, x, y):
    u = coprime_part(gcd(x,y), p*d)
    assert u < 3 or u**3*gcd(d,x*y)**2 < 27*d*d, (p,d,x,y,u)


def check_splits(c, p, gap):
    k = p**gap
    total = (k+1,c)
    q = p
    finite = 0
    while q <= c:
        d = (k+1)*q+c
        for x in range(1,d-1):
            for y in range(x,d-x):
                if low_digits_fit(x,y,d-x-y,p):
                    exclude_numeric(p,d,x,y)
                    finite += 1
        q *= p
    first_q = q
    rays = set()
    for x, y, z in split_forms(c,p,k):
        dx, dy = abs(determinant(total,x)), abs(determinant(total,y))
        if not dx or not dy:
            proportional = x if not dx else y
            assert coprime_part(proportional[0],p) < 3
            continue
        if x[0] and y[0] and determinant(x,y) == 0:
            s, t = y[1], z[1]
            assert x == (k,k*s) and y == (1,s) and z == (0,t)
            assert s > 0 and t > 0 and gcd(p,c) == 1
            assert coprime_part(s,p) == 1
            b = 0
            while p**b != s:
                b += 1
            if coprime_part(t,p) == 1:
                v = 0
                while p**v != t:
                    v += 1
                rays.add((b,v))
                continue
            assert s == 1
            q = first_q
            while (27*(k+1)**2)**2*q <= t**4*(t+1)**6:
                exclude_numeric(p, (k+1)*q+c, k*(q+s), q+s)
                finite += 1
                q *= p
            continue
        u_bound = abs(determinant(x,y)) if x[0] or y[0] else gcd(x[1],y[1])
        assert u_bound > 0
        q = first_q
        while 27*((k+1)*q+c)**2 <= u_bound**3*(dx*dy)**2:
            exclude_numeric(p, (k+1)*q+c, x[0]*q+x[1], y[0]*q+y[1])
            finite += 1
            q *= p
    assert rays == RAYS.get((c,p,gap),set())
    return finite


def pair_carries(x, y, p):
    while x and y:
        if x%p + y%p >= p:
            return True
        x //= p
        y //= p
    return False


def obstruction(c, p, k, b, v, d, g, row, atoms):
    must = g//gcd(*row)
    if gcd(must,k+1) > 1:
        return 'unshifted pair'
    required = set(atoms)
    if p != 2:
        required.add(2)
    for r in required:
        for i in range(3):
            x, y = row[i-1], row[i]
            if (x+y)//gcd(x,y) % r == 0:
                return 'companion pair'
    for r, w in atoms.items():
        tails = [z//w * pow(g,-1,w) % w for z in row if z%w == 0]
        assert len(tails) == 2
        if pair_carries(*tails,r):
            return 'companion digits'
    for r in factors(must):
        unit = coprime_part(g,r)
        w = g//unit
        digit = -p**b * unit * pow(p**v,-1,w) % w
        if pair_carries(digit,k*digit%w,r):
            return 'first-row digits'
    raise ArithmeticError(('remaining row', c,p,k,b,v,d,g,row))


def check_companions(c,p,gap,b,v):
    counts = Counter()
    k = p**gap
    h = 0
    while p**(h+1) <= 2**BIT_BOUND:
        h += 1
        q = p**h
        if q <= c:
            continue
        g = p**(h-b)+1
        d = (k+1)*q+c
        if 27*d*d > g**3:
            continue
        atoms = factors(d)
        if 2*max(atoms.values())+1 > g:
            continue
        for row in complement_by_parts(g,tuple(atoms.values())):
            counts[obstruction(c,p,k,b,v,d,g,row,atoms)] += 1
    return counts


if __name__ == '__main__':
    total = Counter()
    for c,p,gap in CASES:
        finite = check_splits(c,p,gap)
        print(f'Offset {c}, base {p}, gap {gap}: {finite} finite allocation checks.',flush=True)
        for b,v in sorted(RAYS.get((c,p,gap),())):
            counts = check_companions(c,p,gap,b,v)
            print(f'  b={b}, v={v}: {sum(counts.values())} rows, {dict(counts)}.',flush=True)
            total.update(counts)
    print(f'Excluded {sum(total.values())} companion rows: {dict(total)}.')
