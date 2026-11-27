"""Verify the leading-digit exclusions by growing the total before the parts."""

import argparse
from math import comb, gcd
from pathlib import Path

from check_five_leading import leading_depth


def no_prefix(p, caps, denominator, length):
    # Given a total n, each part has an integer maximum. Distribute the small
    # difference from their sum, then test the preceding prefix and digit sum.
    previous = {(0,0,0)}
    for _ in range(length):
        totals = {p*sum(v)+digit for v in previous for digit in range(p)}
        new = set()
        for n in totals-{0}:
            a,b,c = [((n+1)*bound-1)//denominator for bound in caps]
            free = a+b+c-n
            for i in range(min(a,free)+1):
                for j in range(min(b,free-i)+1):
                    x,y,z = a-i,b-j,c-(free-i-j)
                    if z>=0 and x%p+y%p+z%p<p and (x//p,y//p,z//p) in previous:
                        new.add((x,y,z))
        previous = new
    return not previous


def certificates(rows):
    covered = 0
    for k,c,q in rows:
        P = 1 << k
        if (k,c) == (44,32985432719355):
            continue
        if (k,c) == (44,8136949755):
            q = 3**58
        caps = (3*P*q,c*q,(2*P-c)*q+1)
        denom = 5*P*q
        for p in (2,3,5):
            depth = leading_depth(p,caps,denom,max_width=5000 if q==3**58 else 500)
            if depth:
                assert no_prefix(p,caps,denom,depth)
                covered += 1
                break
        else:
            raise AssertionError((k,c))
    assert covered == 31
    # Directly use the fixed low residue at 44 and test all 92 bits by
    # successive factorial quotients, with high zeros at the end.
    k,c,q = 44,8136949755,28147172193075
    P = 1 << k
    valid = {q}
    for j in range(k+2,(3**58).bit_length()+1):
        m = 1 << j
        choices = valid | {v+m//2 for v in valid}
        valid = {v for v in choices if v < 3**58 and 3*v%m+2*v%m<m
                 and (5*P*v+1)//m == (3*P*v)//m+c*v//m+((2*P-c)*v+1)//m}
    assert len(valid) == 7
    for v in valid:
        n,a,b,z = 1+5*P*v,3*P*v,c*v,(2*P-c)*v+1
        assert sum(bin(x).count("1") for x in (a,b,z)) > bin(n).count("1")
    print("Confirmed 31 leading-prefix exclusions and the short binary check at 44.")


def small_controls():
    for p in (2,3,5):
        assert not no_prefix(p,(1,1,1),1,3)
        for n in range(3,48):
            for a in range(1,n-1):
                for b in range(1,n-a):
                    c = n-a-b
                    coefficient = comb(n,a)*comb(n-a,b)
                    if coefficient%p:
                        # Capping each share at its actual value must allow
                        # the three most significant digits, including padding.
                        assert not no_prefix(p,(a,b,c),n,3)
    for d in range(3,102,2):
        for h in range(1,d,2):
            if gcd(h,d)>1:
                continue
            states = range(h-d,h)
            graph = {s:[(bit,b,t) for bit,i,b in [(0,0,0),(1,1,0),(1,0,1),(1,0,0)]
                         if (t:=2*s+d*i-h*bit) in states] for s in states}
            incoming = {s:0 for s in states}
            for edges in graph.values():
                for _,_,t in edges:
                    incoming[t] += 1
            assert max(incoming.values()) <= 2
            def fails(column, word):
                active = set(states)
                for digit in word:
                    active = {t for s in active for x,b,t in graph[s] if (x,b)[column]==digit}
                return not active
            zeros = [0]*d.bit_length()
            assert fails(0,zeros+[1]+zeros) and fails(1,[1]*len(zeros))
    for m in range(1,8):
        P = 1 << m
        for b in range(1,2*P-1):
            n,a,c = 5*P-1,3*P,2*P-1-b
            assert (comb(n,a)*comb(n-a,b))%8 == 4
    print("Exact positive controls and the odd rational slice checks passed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("rows",type=Path,nargs="*")
    args = parser.parse_args()
    small_controls()
    if args.rows:
        certificates([tuple(map(int,line.split())) for f in args.rows for line in f.read_text().splitlines()])
