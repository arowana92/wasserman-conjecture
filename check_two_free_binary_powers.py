"""Every 1+2^a(2^n+2^h+c) for c=83 or 2579, n>h>=1, a>=1."""
from collections import Counter
from itertools import product
from math import gcd, prod

from check_block_ratio_families import possible_companion
from check_endpoint_pair_bound import prime_powers


def certify_shape(n, total, x, y, cap, assigned, counts):
    residue = 1
    caps = set()
    for period in range(1, 4*n+65):
        residue = 2*residue % cap
        u = gcd(total*residue+1,cap)
        if u**3*assigned**2 >= 27*total**2:
            caps.add(u)
        if residue==1:
            counts["closed orbit"] += 1
            break
    else:
        caps = {cap}
    if not caps:
        counts["orbit exclusion"] += 1
        return
    counts["companion stage"] += 1
    forced = [w for p,w in prime_powers(total) if not (
        (x%w==0 and (y%w==0 or y%p!=0)) or (y%w==0 and x%p!=0))]
    for u in caps:
        if 27*prod(forced)**2 <= u**3 and all(2*w<u for w in forced):
            assert not possible_companion(u,forced), (n,total,x,y,u,forced)


def search(c):
    counts = Counter()
    occupied = [i for i in range(c.bit_length()) if c>>i & 1]
    differences = [j-i for i in occupied for j in occupied if i<j]
    assert len(occupied)>=4 and len(set(differences))==len(differences)
    constants = []
    for slots in product(range(3),repeat=len(occupied)):
        constants.append(tuple(sum(1<<i for i,k in zip(occupied,slots) if k==j)
                               for j in range(2)))
    max_n = ((2*c)**64).bit_length()-1
    for n in range(2,max_n+1):
        high = 2**n
        for h in range(1,n):
            low = 2**h
            D = high+low+c
            target = 27*D*D
            if target >= (low+c)**5:
                continue
            B = high//low+1
            if low>c and target>(B*c)**7:
                match = [i for i in occupied if i+n-h in occupied]
                if not match:
                    continue
                R = c-(B<<match[0])
                if target**2 > low**3*(R+1)**6*R**4:
                    continue
            counts["weights"] += 1
            for fx,fy in ((high+low,0),(0,high+low),(0,0),(high,low),(low,high),
                          (high,0),(low,0),(0,high),(0,low)):
                for i,j in constants:
                    x,y = fx+i, fy+j
                    if not 0<x<=y or x+y>=D:
                        continue
                    counts["assignments"] += 1
                    cap = gcd(x,y)
                    cap //= cap & -cap
                    if cap<3:
                        continue
                    assigned = gcd(D,x*y)
                    if cap**3*assigned**2<target:
                        continue
                    g = gcd(cap,D)
                    while g>1:
                        cap //= g
                        g = gcd(cap,D)
                    if cap**3*assigned**2<target or (x&y) or ((x|y)&(D-x-y)):
                        continue
                    counts["candidates"] += 1
                    certify_shape(n,D,x,y,cap,assigned,counts)
    return max_n,counts


if __name__=="__main__":
    for c,expected in ((83,(472,37349,10084279,114,114,114,0)),
                       (2579,(789,105770,95510427,511,505,501,10))):
        cutoff,counts = search(c)
        result = (cutoff, *(counts[k] for k in ("weights","assignments","candidates",
                                                "closed orbit","orbit exclusion","companion stage")))
        print(c,result,flush=True)
        assert result==expected
