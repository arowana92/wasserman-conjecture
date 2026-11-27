"""Test all lifts at an assigned prime, and cover bounded exponent progressions."""

from argparse import ArgumentParser
from collections import Counter, deque
from pathlib import Path


def period_depth(base, step):
    mod = base
    while pow(2,step,mod*base)==1:
        mod *= base
    return mod


def low_ok(slopes, P, modulus, base):
    x,y,z = slopes
    while modulus>1:
        if x*P%modulus + y*P%modulus + (z*P+1)%modulus >= modulus:
            return False
        modulus //= base
    return True


def graph(slopes, base, exponent, step, assigned, cofactor, cap):
    mod = period_depth(base,step)
    # At an assigned weight prime insist also on the exact valuation.
    while assigned and mod<=assigned:
        mod *= base
    initial = pow(2,exponent,mod)
    factor = pow(2,step,mod)
    states = set()
    P = initial
    x,y,z = slopes
    while True:
        if ((not assigned or (cofactor*P+1)%(base*assigned)) and low_ok(slopes,P,mod,base)):
            states.add((x*P//mod,y*P//mod,(z*P+1)//mod))
        P = P*factor%mod
        if P == initial:
            break
    if not states:
        return False
    if base>cap:
        return True
    pending = deque(states)
    while pending:
        c1,c2,c3 = pending.popleft()
        for d in range(base):
            u,v,w = x*d+c1,y*d+c2,z*d+c3
            if u%base+v%base+w%base>=base:
                continue
            nxt = (u//base,v//base,w//base)
            if nxt==(0,0,0) or len(states)>=cap:
                return True
            if nxt not in states:
                pending.append(nxt)
                states.add(nxt)
    return False


def sieve_line(row, bound, counts, budget=10000):
    p,q,a,b,r,w,h,step = row
    slopes = (p*a,q*b,p*q-p*a-q*b)
    n0 = h or step
    if n0>bound:
        return []
    # Retain residue classes modulo increasing powers of r. After the step of
    # a class exceeds the remaining interval, its least exponent is the sole candidate.
    modulus = period_depth(r,step)
    todo = [(n0,step)]
    out = []
    while todo:
        later = []
        for n,H in todo:
            if n+H>bound:
                out.append(n)
                continue
            newmod = modulus*r
            P = pow(2,n,newmod)
            f = pow(2,H,newmod)
            for digit in range(r):
                e = n+H*digit
                if e>bound:
                    break
                if low_ok(slopes,P,newmod,r) and (newmod<=w or (p*q*P+1)%(r*w)):
                    later.append((e,H*r))
                P = P*f%newmod
        modulus *= r
        todo = later
        counts["exponent classes"] += len(todo)
        if len(todo)>budget:
            # Each class can instead be exhausted directly if the range is small.
            if sum((bound-e)//H+1 for e,H in todo)<=budget*4:
                for e,H in todo:
                    out.extend(range(e,bound+1,H))
                break
            raise RuntimeError(("increase covering budget",row,len(todo)))
    return out


def check_candidate(row, k):
    p,q,a,b,r,w,_,_ = row
    slopes = p*a,q*b,p*q-p*a-q*b
    # Doubling the number of layers reaches the actual coefficient if necessary.
    nbits = k+(p*q).bit_length()+1
    levels = 4
    while True:
        complete = True
        for base in (r,p,q):
            mod = base**levels
            if mod.bit_length()<=nbits:
                complete = False
            P = pow(2,k,mod)
            if not low_ok(slopes,P,mod,base):
                return
            if base==r and mod>w and (p*q*P+1)%(r*w)==0:
                return
        if complete:
            raise ArithmeticError(("check the companion",row,k))
        levels *= 2


if __name__=="__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("lines",type=Path)
    parser.add_argument("--cap",type=int,default=1000)
    parser.add_argument("--bound",type=int,default=10**12)
    parser.add_argument("--exclusions", type=Path, help="write the rejecting base for each excluded line")
    args = parser.parse_args()
    counts = Counter()
    surviving = []
    excluded = []
    for line in args.lines.read_text().splitlines():
        if not line or not line[0].isdigit():
            continue
        row = tuple(map(int,line.split()))
        p,q,a,b,r,w,e,H = row
        slopes=(p*a,q*b,p*q-p*a-q*b)
        for base in (r,p,q):
            if not graph(slopes,base,e,H,w if r==base else 0,p*q,args.cap):
                excluded.append((*row,base))
                counts["excluded by base "+("weight" if base==r else ("p" if base==p else "q"))] += 1
                break
        else:
            surviving.append(row)
            if args.bound:
                for k in sieve_line(row,args.bound,counts):
                    counts["bounded exponents checked"] += 1
                    check_candidate(row,k)
    if args.exclusions:
        args.exclusions.write_text("".join(" ".join(map(str,entry))+"\n" for entry in excluded))
    print(dict(counts))
    print("remaining pairs",len({s[:2] for s in surviving}),"progressions",len(surviving))
    for row in surviving:
        print(*row)
