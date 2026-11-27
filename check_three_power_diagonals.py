"""Enumerate close exponent differences by low digit splits and by divisors."""

from argparse import ArgumentParser
from collections import Counter
from itertools import combinations, product
from math import gcd, prod

from search_endpoint import smallest_factors


def shape_works(l, h, k, shifted):
    row = (l, l << h, l << k)
    slots = range(3) if shifted else (-1,)
    for slot in slots:
        a, b, c = [v + (i == slot) for i, v in enumerate(row)]
        if (a | b | c) == a+b+c:
            return True
    return False


def by_scales(limit):
    counts = bytearray(limit+1)
    for k in range(2, limit.bit_length()):
        for h in range(1,k):
            s = 1 + (1<<h) + (1<<k)
            for l in range(1,limit//s+1):
                c = s*l
                if shape_works(l,h,k,False):
                    counts[c] += 1
                if c < limit and shape_works(l,h,k,True):
                    counts[c+1] += 1
    return counts


def by_divisors(counts):
    limit = len(counts)-1
    table = smallest_factors(limit)
    def divisors(n):
        ds = [1]
        while n > 1:
            p = table[n]
            old = ds[:]
            power = 1
            while n % p == 0:
                n //= p
                power *= p
                ds.extend(d*power for d in old)
        return ds
    for c in range(2,limit+1):
        count = 0
        for shifted in (0,1):
            for s in divisors(c-shifted):
                if s & 1 == 0:
                    continue
                rest = s-1
                hbit = rest & -rest
                kbit = rest - hbit
                if hbit < 2 or kbit <= hbit or kbit & (kbit-1):
                    continue
                l = (c-shifted)//s
                parts = (l,l*hbit,l*kbit)
                for slot in (range(3) if shifted else (-1,)):
                    a,b,d = [x+(i==slot) for i,x in enumerate(parts)]
                    if not (a & b or a & d or b & d):
                        count += 1
                        break
        assert counts[c] == count, c


def check_common_factor():
    checks = Counter()
    for p, t, g in product((2,3,5),(1,3,5,7),(3,5,7,11)):
        if t % p == 0:
            continue
        z = p**g
        for indices in ((1,1,1),(3,3,1),(5,3,1),(5,3,2),(4,3,2),(7,5,3)):
            for signs in product((-1,1), repeat=3):
                w = gcd(*(z**i-sign*t**i for i,sign in zip(indices,signs)))
                if w > 2:
                    assert (z-t) % w == 0 or (z+t) % w == 0
                    if indices == (1,1,1):
                        assert len(set(signs)) == 1
                    checks["growing divisor" if w >= z-7 else "other divisor"] += 1
    # After fixing the rational exponent ratios, at most one bin has the
    # required positive absolute values. Check all four targets independently.
    for row in product(range(1,14),repeat=3):
        for degrees in ((1,1,1),(3,2,1),(2,1,1)):
            for shifts in ((0,0,0),(1,0,2)):
                bins = 0
                for e in (-1,0,1,2):
                    x=[abs(v-int(i==e)) for i,v in enumerate(row)]
                    if not all(x):
                        continue
                    first = x[0]
                    ok = all(first**degrees[i] * 2**(shifts[i]*degrees[0]) ==
                             x[i]**degrees[0] * 2**(shifts[0]*degrees[i]) for i in (1,2))
                    bins += ok
                assert bins <= 1
                checks["unique target"] += 1
    return checks


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--limit", type=int, default=2**20)
    args = parser.parse_args()
    data = by_scales(args.limit)
    by_divisors(data)
    print("Algebra checks:", dict(check_common_factor()))
    for limit in sorted({v for v in (256,4096,65536,args.limit) if v<=args.limit}):
        freq = Counter(data[2:limit+1])
        print(limit, "offsets excluded:", freq[0], "remaining:", limit-1-freq[0],
              "offset/difference pairs:", sum(k*n for k,n in freq.items()),
              "maximum differences:", max(freq))
