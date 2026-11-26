"""Recover the companion rows by prime-power partitions and test stable residue layers."""

from itertools import combinations
from math import gcd

from check_close_power_gaps import CASES, RAYS, BIT_BOUND
from check_six_bit_rays import normalized_rows
from check_six_digit_gaps import factors
from check_three_factors import valuation


def full_layer(n, r):
    w = 1
    while n % (w*r) == 0:
        w *= r
    return w


def cross_companions(c, p, gap, b, v):
    k = p**gap
    count = 0
    for h in range(b+1,BIT_BOUND+1):
        q = p**h
        if q <= c or q > 2**BIT_BOUND:
            continue
        x = p**(h-b)
        g = x+1
        d = (k+1)*q+c
        if 27*d*d > g**3:
            continue
        primes = factors(d)
        if max(primes.values()) >= (g+1)//2:
            # For even g the power g/2 is already incompatible with a third positive part.
            continue
        # Evaluate at the first positive outer exponent of the required parity class.
        odd = 1
        while odd*(h-b) <= v:
            odd += 2
        scale = p**(odd*(h-b)-v)
        n = d*scale+1
        first = (k*(q+p**b)*scale,(q+p**b)*scale,p**v*scale+1)
        assert sum(first) == n and all(t%g == 0 for t in first)
        alternate = scale*x*x
        for row in normalized_rows(g,tuple(primes.values())):
            count += 1
            other = tuple(n//g*t for t in row)
            must = g//gcd(*row)
            # A carry between two parts of the same r-adic order, with their
            # unit residues opposite, is permanent under every common scaling.
            forced = set(primes)
            if p != 2:
                forced.add(2)
            if any((z+w)//gcd(z,w)%r == 0 for z,w in combinations(row,2) for r in forced):
                continue
            if any((k+1)%r == 0 for r in factors(must)):
                continue
            failed = False
            for r,w in primes.items():
                indices = [i for i in range(3) if row[i]%w == 0]
                u1,u2 = [other[i] for i in indices]
                assert all(((d*alternate+1)//g*row[i]-other[i])%(w*w) == 0 for i in indices)
                layer = w*r
                while layer <= w*w:
                    failed |= u1%layer + u2%layer >= layer
                    layer *= r
            for r in factors(must):
                w = full_layer(g,r)
                assert all(z*(x*x-1)%(w*w) == 0 for z in first[:2])
                layer = w*r
                while layer <= w*w:
                    failed |= first[0]%layer + first[1]%layer >= layer
                    layer *= r
            assert failed, (c,p,gap,b,v,h,row)
    return count


def arithmetic_controls():
    full = moving = 0
    for c,p,gap in CASES:
        k = p**gap
        for h in range(1,21):
            q = p**h
            d = (k+1)*q+c
            # Enumerate carry-free triples by assigning actual digit units.
            allocations = {(0,0)}
            place, rest = 1,d
            while rest:
                rest, digit = divmod(rest,p)
                allocations = {(x+i*place,y+j*place) for x,y in allocations
                               for i in range(digit+1) for j in range(digit-i+1)}
                place *= p
            for x,y in allocations:
                z = d-x-y
                if not 0 < y <= x or z < 1:
                    continue
                assert valuation(d,(x,y,z),p) == 0
                full += 1
                common = gcd(x,y)
                while gcd(common,p*d) > 1:
                    common //= gcd(common,p*d)
                if common < 3 or common**3*gcd(d,x*y)**2 < 27*d*d:
                    continue
                # Use the actual weight across a full parity period for each
                # proportional row, and independently apply the size inequality.
                s = x//k-q
                assert x == k*y and y == q+s and s > 0
                b = 0
                while p**b < s:
                    b += 1
                assert p**b == s and h > b
                for a in range(1,2*(h-b)+1):
                    u = gcd(common,1+p**a*d)
                    if u**3*gcd(d,x*y)**2 >= 27*d*d:
                        options = RAYS.get((c,p,gap),set())
                        assert any(i == b and (a+v)%(2*(h-b)) == h-b for i,v in options)
                        moving += 1
    print(f'{full} finite digit allocations and {moving} exponent controls.',flush=True)


if __name__ == '__main__':
    arithmetic_controls()
    total = 0
    for (c,p,gap), rays in RAYS.items():
        for b,v in sorted(rays):
            count = cross_companions(c,p,gap,b,v)
            print(f'{(c,p,gap,b,v)}: {count} rows by prime-power partitions.',flush=True)
            total += count
    print(f'Independently excluded {total} companion rows.')
