"""Independent cofactor and modular-gcd check for the two binary offsets."""
from collections import Counter
from math import gcd

from check_block_ratio_families_cross import companion_exists, powers


def check_cofactor(D, n, counts):
    need = 27*D**2
    union = (D-1)&D
    while union:
        free = union^(union & -union)
        cap_union = union//(union & -union)
        common = gcd(cap_union,D)
        while common>1:
            cap_union //= common
            common = gcd(cap_union,D)
        x = free if cap_union>=3 else 0
        while x:
            y = union-x
            cap = gcd(x,cap_union)
            if cap>=3:
                first = gcd(D,x)
                assigned = first*gcd(D//first,y)
                if cap**3*assigned**2 >= need:
                    counts["product stage"] += 1
                    for j in range(1,2*n+65):
                        m = gcd(cap,(1<<j)-1)
                        if m>1:
                            cap = (cap//m)*gcd(m,pow(-D,j,m)-1)
                        if cap**3*assigned**2 < need:
                            break
                    if cap**3*assigned**2 >= need:
                        counts["subsets"] += 1
                        forced = [w for p,w in powers(D) if (x%w and y%w) or
                                  (x%w and x%p==0) or (y%w and y%p==0)]
                        assert not companion_exists(cap,forced), (D,x,y,cap,forced)
            x = (x-1)&free
        union = (union-1)&D


def search(c):
    counts = Counter()
    positions = [i for i in range(c.bit_length()) if c&(1<<i)]
    distances = {v-u:u for u in positions for v in positions if v>u}
    assert 2*len(distances)==len(positions)*(len(positions)-1)
    n = 2
    while 2**n <= (2*c)**64:
        for gap in range(1,n):
            P = 2**(n-gap)
            B = 2**gap+1
            D = P*B+c
            if (P+c)**5 <= 27*D**2:
                continue
            if P>c and (B*c)**7 < 27*D**2:
                if gap not in distances:
                    continue
                rest = c-B*2**distances[gap]
                if P**3*rest**4*(rest+1)**6 < 729*D**4:
                    continue
            counts["weights"] += 1
            check_cofactor(D,n,counts)
        n += 1
    return counts


if __name__=="__main__":
    for c,expected in ((83,(37349,114,3)), (2579,(105770,511,15))):
        counts = search(c)
        print(c,counts,flush=True)
        assert tuple(counts[k] for k in ("weights","product stage","subsets")) == expected
