"""Check predecessor defects directly using the digits of the unmodified row."""

from collections import deque
from itertools import combinations
from math import comb, gcd, lcm

from check_affine_endpoint import factors
from check_direction_predecessor_loss import constant, prime_loss, vp


def raw_accepts(row,p,e,j):
    q=p**j
    prefix=pow(sum(row),-1,q)
    residues=[(a*prefix)%q for a in row]
    assert sum(residues)==1
    # The initial suffix of N is 1 followed by zeros. Each further such
    # zero earns a unit of predecessor valuation until the first mismatch.
    start=(tuple(a*prefix//q for a in row),0,True,0)
    todo=deque([start])
    seen={start}
    while todo:
        pending,carry,matching,credit=todo.popleft()
        if pending==(0,0,0) and carry==0 and not matching:
            return True
        for digit in range(p):
            parts=tuple(x+a*digit for x,a in zip(pending,row))
            outgoing,n_digit=divmod(sum(v%p for v in parts)+carry,p)
            continues=matching and n_digit==0
            balance=credit+int(continues)-outgoing
            if balance<0:
                continue
            assert balance<=e-j
            dest=(tuple(v//p for v in parts),outgoing,continues,balance)
            if dest not in seen:
                seen.add(dest)
                todo.append(dest)
    return False


def run(limit):
    raw=exact=reduced=0
    for x in range(1,limit//3+1):
        for y in range(x,(limit-x)//2+1):
            for z in range(y,limit-x-y+1):
                if gcd(x,y,z)!=1:
                    continue
                row=(x,y,z)
                h=lcm(*(gcd(a,b) for a,b in combinations(row,2)))
                c=constant(row)
                assert h%c==0
                if c<h:
                    reduced+=1
                for p,w in factors(h):
                    e=vp(w,p)
                    expected,_=prime_loss(row,p,w)
                    for j in range(1,e+1):
                        assert raw_accepts(row,p,e,j)==(expected>=p**j)
                        raw+=1
                if sum(row)<=25:
                    for v in range(1,26):
                        n=sum(row)*v
                        coefficient=comb(n,x*v)*comb(n-x*v,y*v)
                        assert c*coefficient%(n-1)==0
                        exact+=1
    print("sum limit",limit,"raw prime-power comparisons",raw,"exact coefficients",exact,
          "strict reductions",reduced,flush=True)


if __name__=="__main__":
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--limit",type=int,default=90)
    run(parser.parse_args().limit)
