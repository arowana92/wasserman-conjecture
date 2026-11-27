"""Exact finite controls for the block ratio and companion arguments."""
from fractions import Fraction
from itertools import combinations, product
from math import comb, gcd

from check_block_ratio_families import possible_companion


def splits(n, p):
    result = []
    for x in range(n+1):
        for y in range(n-x+1):
            if comb(n,x)*comb(n-x,y) % p:
                result.append((x,y,n-x-y))
    return result


def minors():
    cases = common_content = 0
    for p in (2,3,5,7):
        tables = {n:splits(n,p) for n in range(1,17)}
        ratios = {n:{Fraction(x,y) for x,y,z in v if x and y} for n,v in tables.items()}
        for b,c in product(tables,repeat=2):
            if not ratios[b].isdisjoint(ratios[c]):
                continue
            P = p
            while P<=c:
                P *= p
            D = b*P+c
            for upper,lower in product(tables[b],tables[c]):
                parts = [i*P+j for i,j in zip(upper,lower)]
                if min(parts)==0:
                    continue
                for i,j,x in zip(upper,lower,parts):
                    assert b*j != c*i and gcd(x,D)<=b*c
                for x,y in combinations(parts,2):
                    cap = gcd(x,y)
                    while cap%p==0:
                        cap //= p
                    assert cap<=b*c
                    assert cap**3*gcd(D,x*y)**2 <= (b*c)**7
                cases += 1
                common_content += gcd(b,c)>1
    assert common_content
    print("Exact carry-free block splits",cases,"with noncoprime block values",common_content)


def matched_direction():
    checks = nontrivial = 0
    for p in (2,3,5,7):
        for v in range(3):
            positions = [v,v+1,v+4,v+6]
            c = sum(p**j for j in positions)
            n = v+8
            P = p**n
            B = p+1
            total = B*P+c
            R = c-B*p**v
            for allocation in product(range(3), repeat=6):
                upper = [0,0,0]
                lower = [0,0,0]
                for digit,slot in zip((1,p),allocation[:2]):
                    upper[slot] += digit
                for pos,slot in zip(positions,allocation[2:]):
                    lower[slot] += p**pos
                if min(*upper[:2],*lower[:2])==0 or upper[0]*lower[1]!=upper[1]*lower[0]:
                    continue
                x,y = sorted(upper[i]*P+lower[i] for i in (0,1))
                assert x==P+p**v and y==p*x
                for a in range(1,80):
                    u = gcd(1+p**a*total,gcd(x,y))
                    assert (p**(n-v)+1)%u==0 and (p**a*R+1)%u==0
                    assert u*u <= P*(R+1)**2
                    checks += 1
                    nontrivial += u>1
    assert nontrivial
    print("Single shared direction checks",checks,"nontrivial scale gcds",nontrivial)


def companion_cases():
    yes = no = 0
    for atoms in ([],[3],[4,5],[3,5,7],[4,7,9]):
        for cap in range(3,65):
            if any(gcd(cap,w)>1 for w in atoms):
                continue
            actual = any(all(sorted([x%w,y%w,(cap-x-y)%w])==[0,0,cap%w]
                             for w in atoms)
                         for x in range(1,cap//3+1) for y in range(x,(cap-x)//2+1))
            assert possible_companion(cap,atoms)==actual
            yes += actual
            no += not actual
    print("Companion controls",yes,"feasible,",no,"impossible")


def quantitative_gcd():
    checked = 0
    for p in (2,3,5):
        for b,c in ((1,1),(1,2),(2,1),(3,1),(1,3)):
            for h in (1,2,3):
                P = p**h
                if P<=c:
                    continue
                D = b*P+c
                for a in (1,2):
                    power = p**a
                    for shift in (0,1):
                        n = power*D+shift
                        if n>125 or (shift==0 and c%p==0):
                            continue
                        coeffs = [comb(n,x)*comb(n-x,y) for x in range(1,n//3+1)
                                  for y in range(x,(n-x)//2+1)]
                        left = {gcd(n*(n-1),k) for k in coeffs if k%p}
                        right = {gcd(n*(n-1),k) for k in coeffs}
                        for first,second in product(left,right):
                            common = gcd(first,second)
                            if shift:
                                assert common**3*(b*c)**7 >= 27*D**2
                            else:
                                assert common**3*power**3*(b*c)**6 >= 27*(n-1)**2
                            checked += 1
    assert checked>0
    print("Quantitative pair gcd controls",checked)


if __name__=="__main__":
    minors()
    matched_direction()
    companion_cases()
    quantitative_gcd()
