"""Check power congruences and Kummer carries for the binary shapes at 1+2^k*p*q."""

from argparse import ArgumentParser
from collections import Counter, defaultdict
from functools import lru_cache
from math import gcd
from pathlib import Path


def power_logs(mod):
    logs = [-1]*mod
    x, step = 1, 0
    while logs[x] == -1:
        logs[x] = step
        x = 2*x % mod
        step += 1
    return logs, step


small_logs = lru_cache(None)(power_logs)


def join(a, m, b, n):
    common = gcd(m,n)
    if (a-b)%common:
        return None
    factor = n//common
    s = ((b-a)//common)*pow(m//common,-1,factor)%factor
    return a+s*m, m*factor


def odd_prime_powers(n):
    while n%2==0:
        n //= 2
    r = 3
    while r*r<=n:
        power=1
        while n%r==0:
            n //= r
            power *= r
            yield r,power
        r += 2
    if n>1:
        yield n,n


def layer_ok(x,y,z,w,Q):
    return (Q*x%w + Q*y%w + (Q*z+1)%w) < w


def finite_ok(x,y,z,p,k):
    scale = 1<<k
    row = (scale*x,scale*y,scale*z+1)
    total = sum(row)
    power = p
    while power <= total:
        if sum(a%power for a in row)>=power:
            return False
        power *= p
    return True


def scan(records, finite):
    counts = Counter()
    possible = set()
    by_q = defaultdict(list)
    for p,q,a,b in records:
        by_q[q].append((p,a,b))
    for q,rows in by_q.items():
        qlog,qo = power_logs(q)
        for p,a,b in rows:
            m,x,y = p*q,p*a,q*b
            z = m-x-y
            assert z>0 and (x|y|z)==m and not (x&y or x&z or y&z)
            plog,po = small_logs(p)
            px,qy = pow(y,-1,p),pow(x,-1,q)
            assert plog[px]>=0
            counts["binary shapes"] += 1
            if qlog[qy]<0:
                continue
            pair = join(plog[px],po,qlog[qy],qo)
            if pair is None:
                continue
            counts["compatible p and q powers"] += 1
            if px*a%p > px*q%p or qy*b%q > qy*p%q:
                continue
            counts["second digits at p and q"] += 1
            for r,w in odd_prime_powers(gcd(a,b)):
                logs,order = small_logs(w)
                required = logs[-pow(m,-1,w)%w]
                if required<0:
                    continue
                line = join(*pair,required,order)
                if line is None:
                    continue
                u,H = line
                R = -pow(m,-1,r)%r
                if (x//w*R)%r + (y//w*R)%r >= r:
                    continue
                counts["weight prime progressions"] += 1
                possible.add((p,q,a,b,r,w,u,H))
                for k in range(u or H,finite+1,H):
                    counts["finite exponents"] += 1
                    if (finite_ok(x,y,z,r,k) and finite_ok(x,y,z,p,k) and finite_ok(x,y,z,q,k)):
                        raise ArithmeticError(("check the companion",p,q,a,b,r,w,k))
    return counts,possible


if __name__=="__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("shapes",type=Path)
    parser.add_argument("--exponents",type=int,default=4096)
    args = parser.parse_args()
    records = [tuple(map(int,line.split())) for line in args.shapes.read_text().splitlines()]
    counts,possible = scan(records,args.exponents)
    print(dict(counts))
    print("remaining pairs",len({s[:2] for s in possible}),"progressions",len(possible))
    for s in sorted(possible):
        print(*s)
