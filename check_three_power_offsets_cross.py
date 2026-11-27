"""Independently enumerate finite candidate rows via a single large target divisor."""

from argparse import ArgumentParser
from collections import Counter
from functools import lru_cache
from math import gcd, isqrt, prod


def retained(n, divisor):
    other = n // divisor
    while (shared := gcd(divisor,other)) > 1:
        other *= shared
        divisor //= shared
    return divisor


def unit_shape(row, p):
    base = min(row)
    for x in row:
        ratio, rest = divmod(x, base)
        if rest:
            return False
        while ratio % p == 0:
            ratio //= p
        if ratio != 1:
            return False
    return True


def digit_sum(n, p):
    total = 0
    while n:
        n, d = divmod(n,p)
        total += d
    return total


def cross(offset, p):
    assert 2 <= offset <= 2048 and offset*(offset-1) % p == 0
    assert (2*offset)**5 < 27*10**18
    height = 9
    scale = 1
    while scale < 2*offset+2:
        height += 90
        scale *= p
    start = 0
    while (p**start + offset)**5 <= 27*10**18:
        start += 1
    assert p**start > offset
    power = [p**i for i in range(height+1)]

    ppart = 1
    while offset*(offset-1) % (p*ppart) == 0:
        ppart *= p
    digits = [digit_sum(x,p) for x in range(offset+1)]
    allocations = [(i,j,offset-i-j) for i in range(offset+1) for j in range(offset-i+1)
                   if digits[i]+digits[j]+digits[offset-i-j] == digits[offset]]
    counts = Counter()
    remaining = set()
    for a in range(start, height+1):
        # One target has ninth power at least 27(N-1)^2, and its p-part is at most ppart.
        target = 27*max(10**9,power[a])**2
        lo, hi = 0, 1 << ((target.bit_length()+8)//9)
        while lo+1 < hi:
            mid = (lo+hi)//2
            if (ppart*mid)**9 >= target:
                hi = mid
            else:
                lo = mid
        threshold = hi

        @lru_cache(None)
        def choices(high_offset, low_offset):
            large = power[a]+high_offset
            while large % p == 0:
                large //= p
            if large < threshold:
                return ()
            return tuple((b,g) for b in range(start,a+int(p!=2))
                         if (g := gcd(large,power[b]+low_offset)) >= threshold)

        visited = set()
        for lowrow in allocations:
            for slot in (-1,0,1,2):
                i,j,k = (v-int(index==slot) for index,v in enumerate(lowrow))
                for b, gb in choices(i,j):
                    for c, gc in choices(i,k):
                        if c>b or (c==b and (p==2 or (p==3 and a==b))):
                            continue
                        if gcd(gb,gc) < threshold or (lowrow,b,c) in visited:
                            continue
                        visited.add((lowrow,b,c))
                        n = power[a]+power[b]+power[c]+offset
                        if n <= 10**9:
                            continue
                        counts["candidate rows"] += 1
                        row = [power[e]+x for e,x in zip((a,b,c),lowrow)]
                        if unit_shape(row,p):
                            counts["unshifted shapes"] += 1
                            continue
                        caps = [gcd(*row)]
                        for j in range(3):
                            adjusted = [x-int(index==j) for index,x in enumerate(row)]
                            if unit_shape(adjusted,p):
                                counts["shifted shapes"] += 1
                                break
                            caps.append(gcd(*adjusted))
                        else:
                            u = retained(n,caps[0])
                            s = retained(n-1,prod(caps[1:]))
                            if u>=3 and n//u>=3 and u**3*s**2 >= 27*(n-1)**2:
                                remaining.add((n,tuple(row)))
        choices.cache_clear()
    return counts, remaining


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("offsets", nargs="+", type=int)
    args = parser.parse_args()
    for offset in args.offsets:
        for p in range(2,offset+1):
            if offset*(offset-1) % p or any(p%d == 0 for d in range(2,isqrt(p)+1)):
                continue
            counts, remaining = cross(offset,p)
            print(offset,p,dict(counts),"unresolved:",len(remaining),flush=True)
            if remaining:
                print(sorted(remaining))
