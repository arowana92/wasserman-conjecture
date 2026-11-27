"""Check exponent separation directly with factorial floors and sorted residue lists."""
from argparse import ArgumentParser
from itertools import combinations_with_replacement
from math import gcd
from typing import NamedTuple


class Row(NamedTuple):
    slopes: tuple
    first: int
    stride: int

    @property
    def d(self):
        return sum(self.slopes)


def factor_five(q, z, first, stride):
    return Row((3*q,z,2*q-z),first,stride)


LARGE = (factor_five(1467779,2918400,2337019,2935556),
         factor_five(2823523,5636100,3948614,5647044))
SMALL = (Row((681,9222,2128),769251571,2043279108),
         Row((16536,36930,11777),174688,287820))


def keep_neighbors(rows, others, period, gap):
    # Scan all entries once. The previous and next *strictly* different
    # residues of the other list are the only ones needed, including wrap.
    if not others:
        return []
    left = right = 0
    kept = []
    for entry in rows:
        e = entry[0]
        while left < len(others) and others[left][0] < e:
            left += 1
        while right < len(others) and others[right][0] <= e:
            right += 1
        prev_e = others[left-1][0] if left else others[-1][0]-period
        next_e = others[right][0] if right < len(others) else others[0][0]+period
        if min(e-prev_e,next_e-e) <= gap:
            kept.append(entry)
    return kept


def next_layer(row, values, period, width, precision, multipliers):
    """Test each possible new exponent digit against the exact row residues."""
    constraint = gcd(3*period,row.stride)
    target = row.first % constraint
    x,y,z = row.slopes
    out = []
    for i, multiplier in enumerate(multipliers):
        for e,q in values:
            n = e+i*period
            if n % constraint != target:
                continue
            value = q*multiplier % precision if i else q
            v = value % width
            if (row.d*v+1) % 9 == 0 and width >= 9:
                continue
            if (x*v)%width+(y*v)%width+(z*v+1)%width < width:
                out.append((n,value))
    return out


def certify(first, second, bound):
    period = gcd(first.stride,second.stride)
    assert period % 2 == 0
    width = 3
    while pow(2,period,3*width) == 1:
        width *= 3
    precision = width*3**360

    def root(row):
        e = row.first % period
        return [(e,pow(2,e,precision))]

    same = first == second
    a, b = root(first), root(second)
    peak = 0
    for depth in range(360):
        if period > bound:
            if same:
                a = keep_neighbors(a,a,period,bound)
            else:
                a,b = keep_neighbors(a,b,period,bound),keep_neighbors(b,a,period,bound)
        if not a or not b:
            return depth,peak
        peak = max(peak,len(a)+(0 if same else len(b)))
        ratio = pow(2,period,precision)
        multipliers = (1,ratio,ratio**2 % precision)
        width *= 3
        a = next_layer(first,a,period,width,precision,multipliers)
        if not same:
            b = next_layer(second,b,period,width,precision,multipliers)
        period *= 3
    raise ArithmeticError(("unseparated",first.d,second.d,bound))


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--mixed-only", action="store_true")
    args = parser.parse_args()
    if not args.mixed_only:
        for r,s in combinations_with_replacement(LARGE,2):
            print(r.d,s.d,128,*certify(r,s,10**128),flush=True)
    for r in LARGE:
        for s in SMALL:
            print(r.d,s.d,15,*certify(r,s,10**15),flush=True)
