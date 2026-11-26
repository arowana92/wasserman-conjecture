"""Certify row-gap pairs with e-d=t-1, de=t*N, at unrestricted weights."""

from argparse import ArgumentParser
from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
from itertools import product
from math import gcd, prod

from check_divisor_cuts import divisors
from search_even_endpoint import prime_powers


def coprime_bins(bins):
    previous = 1
    for value in bins:
        if gcd(value, previous) > 1:
            return False
        previous *= value
    return True


@lru_cache(maxsize=16384)
def multipliers(denominator, middle_bin, last_bin):
    # denominator = 2*m+h, with last_bin | m and middle_bin | h.
    return tuple((m, denominator-2*m)
                 for m in range(last_bin, (denominator-middle_bin)//2+1, last_bin)
                 if (denominator-2*m) % middle_bin == 0
                 and gcd(m, middle_bin) == gcd(denominator-2*m, last_bin) == 1)


@dataclass(frozen=True)
class Shape:
    x: int
    y: int
    b: int
    c: int
    m: int
    h: int
    n: int
    j: int
    s2: int
    s3: int
    t1: int
    t3: int

    def fixed_bins(self):
        return self.s2, self.s3, self.t1, self.t3

    def rows(self, s1, t2):
        u, v = self.c*t2-self.y, self.b*s1+self.x
        a = (u*(self.m*s1+self.x), u*self.m*s1, u*self.h*s1)
        z = (v*self.n*t2, v*(self.n*t2-self.y), v*self.j*t2)
        return u, v, a, z

    def progression(self):
        # x*c*t2 - y*b*s1 = xy+1, starting with the least nonnegative s1.
        left, right = self.x*self.c, self.y*self.b
        common = gcd(left, right)
        s_step, t_step = left//common, right//common
        target = self.x*self.y+1
        s1 = (-(target//common)*pow(t_step, -1, s_step)) % s_step
        t2, rem = divmod(target+right*s1, left)
        assert rem == 0
        return s1, t2, s_step, t_step


def shapes(t):
    for x in divisors(t):
        y = t//x
        if gcd(x, y) > 1:
            continue
        for s3, t3, s2, t1 in product(divisors(x), divisors(y),
                                        divisors(2*(t-1)), divisors(2*(t-1))):
            bins = s2, s3, t1, t3
            if not coprime_bins(bins) or (min(s2,s3) == 1 and min(t1,t3) == 1):
                continue
            size = prod(bins)
            for b in divisors(size):
                c = size//b
                if b < 2*s3+s2 or c < 2*t3+t1 or (t+1) % gcd(x*c, y*b):
                    continue
                for (m,h), (n,j) in product(multipliers(b,s2,s3), multipliers(c,t1,t3)):
                    yield Shape(x,y,b,c,m,h,n,j,s2,s3,t1,t3)


def fails_at_base(shape, p, z, modulus, progression):
    s0, t0, ds, dt = progression
    s1, t2 = s0+ds*z, t0+dt*z
    u, v, a, c = shape.rows(s1,t2)
    bins = (*shape.fixed_bins(), s1, t2)
    if (u % p == v % p == 0) or sum(w % p == 0 for w in bins) > 1:
        return True
    for row, assigned in ((a, u*s1*shape.s2*shape.s3), (c, v*t2*shape.t1*shape.t3)):
        if assigned % p == 0 and sum(part % modulus for part in row) != (u*v) % modulus:
            return True
    return False


def first_nonzero_carry(shape, p, residue, progression):
    """Test the digit at the valuation of a positive variable scale or bin."""
    s0,t0,ds,dt = progression
    s1,t2 = s0+ds*residue,t0+dt*residue
    u,v = shape.c*t2-shape.y,shape.b*s1+shape.x
    m,h,n,j = shape.m,shape.h,shape.n,shape.j
    candidates = (
        (u, (m*s1+shape.x, m*s1, h*s1)),
        (v, (n*t2, n*t2-shape.y, j*t2)),
        (s1, (shape.b*shape.y+u*m, u*m, u*h)),
        (t2, (v*n, v*n-shape.c*shape.x, v*j)),
    )
    for factor, slopes in candidates:
        if factor % p == 0 and all(sum(unit*a % p for a in slopes) >= p
                                   for unit in range(1,p)):
            return True
    return False


def local_witness(shape):
    # Every tested modulus is exact. Exhausting a search cap never proves exclusion.
    progression = shape.progression()
    known = shape.x*shape.y*prod(shape.fixed_bins())
    for p in sorted({2,3} | set(prime_powers(known))):
        blocked = {r for r in range(p) if first_nonzero_carry(shape,p,r,progression)} if p <= 17 else set()
        w = 1
        residues = [0]
        while w*p <= 16384:
            q = w*p
            residues = [r for old in residues for r in range(old,q,w)
                        if r % p not in blocked and not fails_at_base(shape,p,r,q,progression)]
            if not residues:
                return p,q
            if len(residues) > 500:
                break
            w = q
    return None


def verify(limit):
    counts = Counter()
    largest = 0
    witnesses = Counter()
    for t in range(2, limit+1):
        for shape in shapes(t):
            counts['shapes'] += 1
            for first in (True,False):
                other_bins = shape.t1*shape.t3 if first else shape.s2*shape.s3
                if other_bins % 2 == 0:
                    continue
                if not first and shape.j % (2*shape.n) == 0:
                    witness = local_witness(shape)
                    if witness is None:
                        raise ArithmeticError(('unresolved unbounded progression', t, shape))
                    witnesses[witness] += 1
                    continue
                s1,t2,ds,dt = shape.progression()
                while (shape.m*s1 < shape.b*shape.x if first
                       else shape.n*t2 < (shape.c+1)*shape.y):
                    u,v,a,c = shape.rows(s1,t2)
                    if min(s1,t2,*a,*c) > 0:
                        weight = u*v
                        assert sum(a) == sum(c) == weight
                        assert (a[0]-a[1])*(c[0]-c[1]) == t*weight
                        assert (c[0]-c[1])-(a[0]-a[1]) == t-1
                        counts['positive'] += 1
                        largest = max(largest,weight)
                        if gcd(u,v) == 1 and coprime_bins((*shape.fixed_bins(),s1,t2)):
                            counts['coprime'] += 1
                            row = a if first else c
                            if row[0] | row[1] | row[2] == weight:
                                raise ArithmeticError(('requires further carry tests', weight, a, c))
                    s1 += ds
                    t2 += dt
    print('Normalized gap product through',limit,dict(counts))
    print('Largest positive finite weight',largest)
    print('Unbounded progressions by killing (prime, modulus)',sorted(witnesses.items()))


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('limit', type=int, nargs='?', default=2660)
    verify(parser.parse_args().limit)
