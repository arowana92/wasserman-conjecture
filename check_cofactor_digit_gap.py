"""Verify the gcd gap, including small denominators and three- and four-factor rows."""

from collections import Counter
from itertools import combinations, product
from math import comb, gcd, prod


def digits(n, base):
    result = []
    while n:
        n, digit = divmod(n,base)
        result.append(digit)
    return result


def value(coefficients, base):
    result = 0
    for digit in reversed(coefficients):
        result = result*base+digit
    return result


def check_gap(P, D, B, counts):
    ds = digits(D,P)
    c = ds[0]
    if not c:
        return
    s = D//gcd(B,D)
    if s*(c+1)*max(ds) >= P:
        counts["outside small-s range"] += 1
        return
    assert gcd(*ds) == c
    F = D//c
    factor, rem = divmod(B,F)
    y, unit = divmod(factor,P)
    assert rem == 0 and unit == 1 and 0 <= y <= c-2
    assert s == c//gcd(c,1+y*P)
    assert all(b <= (c-y)*a for a,b in zip(ds,ds[1:]))
    assert gcd(B,P*D+1) == gcd(1+y*P,P*D+1)
    counts["small-s conclusions"] += 1
    if s > 1:
        counts["small-s with s > 1"] += 1
    if y and gcd(B,P*D+1) >= 3:
        counts["shared divisor at least 3"] += 1


def row_controls():
    counts = Counter()
    for P in (5,8,9,13,17,29,49,97,257):
        for size in (1,2,3,4):
            for ds in product(range(4),repeat=size):
                if not ds[0] or not ds[-1]:
                    continue
                D = value(ds,P)
                for allocation in product(*(range(x+1) for x in ds)):
                    if sum(ds)-sum(allocation) < 2:
                        continue
                    B = 1+P*value(allocation,P)
                    check_gap(P,D,B,counts)
    # Test the nonconstant-content and nontrivial-denominator branches directly.
    for P in (101,127,157,251,509,977):
        for c in range(3,12):
            for y in range(1,c-1):
                for fs in ((1,1),(1,2,1),(1,1,1,1),(1,3,2)):
                    if any(b > (c-y)*a for a,b in zip(fs,fs[1:])):
                        continue
                    F = value(fs,P)
                    check_gap(P,c*F,(1+y*P)*F,counts)
    assert all(counts[k] for k in ("outside small-s range","small-s conclusions",
                                   "small-s with s > 1","shared divisor at least 3"))
    print("Digit and gcd checks",dict(counts))


def exact_controls():
    count = 0
    for p,P in ((2,8),(2,16),(3,9),(3,27),(5,5),(7,7),(11,11),(17,17)):
        for D in range(1, min(60,450//P)+1):
            n = 1+P*D
            for a in range(1,n//3+1):
                for b in range(a,(n-a)//2+1):
                    row = (a,b,n-a-b)
                    if comb(n,a)*comb(n-a,b) % p == 0:
                        continue
                    special, = [v for v in row if v%P]
                    assert special%P == 1
                    size = len(digits(n,P))
                    assert all(sum(v//P**j%P for v in row) == n//P**j%P
                               for j in range(size))
                    check_gap(P,D,special,Counter())
                    count += 1
    print(f"{count} carry-free rows from exact integer trinomial coefficients.")


def atoms(bound):
    result = []
    for p in range(2,bound+1):
        if any(p%q == 0 for q,_ in result if q*q <= p):
            continue
        w = p
        while w <= bound:
            result.append((p,w))
            w *= p
    return result


def crt_controls():
    counts = Counter()
    for triple in combinations(atoms(95),3):
        if len({p for p,_ in triple}) != 3:
            continue
        for extra in (1,2):
            if extra == 2 and any(p == 2 for p,_ in triple):
                continue
            m = extra*prod(w for _,w in triple)
            for doubled_slot in (range(3) if extra == 2 else (None,)):
                bins = [2*w if i==doubled_slot else w for i,(_,w) in enumerate(triple)]
                # For the four-factor case this gives all placements of 2 on this row.
                parts = [(m//w)*pow(m//w,-1,w) for w in bins]
                if sum(parts) != m+1 or gcd(*parts)<3:
                    continue
                for (p,w),B in zip(triple,parts):
                    P = p
                    while P <= w:
                        D = m//P
                        assert D//gcd(B,D) <= extra*(w//P)
                        if extra == 1:
                            assert D//gcd(B,D) == w//P
                        P *= p
                        counts[extra] += 1
    assert counts[1] and counts[2]
    print(f"{counts[1]} three-factor and {counts[2]} four-factor CRT divisor checks.")


if __name__ == "__main__":
    row_controls()
    exact_controls()
    crt_controls()
