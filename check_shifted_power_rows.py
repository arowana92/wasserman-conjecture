"""Direct coefficient, carry-layer and binary offset checks for shifted rows."""

from argparse import ArgumentParser
from collections import Counter
from functools import lru_cache
from itertools import combinations_with_replacement
from math import comb, gcd

from search_endpoint import smallest_factors


def vp(n,p):
    v=0
    while n%p==0:
        v+=1
        n//=p
    return v


def layers(row,p):
    total=sum(row)
    val=0
    q=p
    while q<=total:
        val+=(sum(z%q for z in row)-total%q)//q
        q*=p
    return val


def carry_checks():
    counts=Counter()
    for p in (2,3,5,7):
        for i in range(1,14):
            for j in range(i,2*i+5):
                s=1+p**i+p**j
                for v in range(1,2*i+6):
                    mod=p**v
                    m=(-pow(s,-1,mod))%mod
                    for step in (0,1):
                        t=m+step*mod
                        row=(t+1,t*p**i,t*p**j)
                        true_v=vp(sum(row),p)
                        lam=layers(row,p)
                        r=max(true_v-lam,0)
                        assert p**(2*r)<=s
                        if p!=2 or i!=j:
                            assert r<=min(i,j-i)
                        assert lam>=max(true_v-i,0)+max(i+min(i,true_v)-j,0)-int(p==2 and i==j)
                        counts["carry profiles"]+=1
                        if 0<r<true_v:
                            counts["partial weight powers"]+=1
                        if j<2*i and true_v>j and i<j:
                            counts["overlapping layers inside weight power"]+=1
    return counts


def exact_pairs(limit=155):
    @lru_cache(None)
    def endpoint_divisors(n):
        return {gcd(n*(n-1),comb(n,a)*comb(n-a,b)) for a in range(1,n//3+1)
                for b in range(a,(n-a)//2+1)}
    counts=Counter()
    for p in (2,3,5,7,11):
        powers=[1]
        while powers[-1]*p<=limit:
            powers.append(powers[-1]*p)
        for y,z in combinations_with_replacement(powers,2):
            direction=(1,y,z)
            s=sum(direction)
            for m in range(1,(limit-1)//s+1):
                n=m*s+1
                for slot in range(3):
                    row=[m*t+int(i==slot) for i,t in enumerate(direction)]
                    a,b,_=row
                    ka=comb(n,a)*comb(n-a,b)
                    for gb in endpoint_divisors(n):
                        g=gcd(ka,gb)
                        assert g>=3
                        if n%p==0:
                            assert g**6 >= 729*s
                        else:
                            primitive=s//p**vp(s,p)
                            assert g**3>=27*primitive**2
                        counts["endpoint gcd comparisons"]+=1
                        if ka%p and n%p==0:
                            counts["designated row avoids weight prime"]+=1
    return counts


def offset_counts(limit):
    candidate=bytearray(limit+1)
    for i in range(1,limit.bit_length()):
        for j in range(i+1,limit.bit_length()):
            s=1+(1<<i)+(1<<j)
            for l in range(1,limit//s+1):
                a,b=l<<i,l<<j
                if not (a&b or a&l or b&l):
                    candidate[l*s]+=1
    spf=smallest_factors(limit)
    def divisors(n):
        ds=[1]
        while n>1:
            p=spf[n]
            layer=ds[:]
            while n%p==0:
                n//=p
                layer=[p*x for x in layer]
                ds.extend(layer)
        return ds
    for c in range(2,limit+1):
        count=0
        for s in divisors(c):
            if s%2==0:
                continue
            p1=(s-1)&-(s-1)
            p2=s-1-p1
            if p1<=1 or p2<=p1 or p2&(p2-1):
                continue
            l=c//s
            if (l|l*p1|l*p2)==c:
                count+=1
        assert candidate[c]==count
    for bound in sorted({b for b in (63,256,4096,65536,limit) if b<=limit}):
        freqs=Counter(candidate[2:bound+1])
        print(bound,"excluded",freqs[0],"remaining",bound-1-freqs[0],
              "pairs",sum(k*v for k,v in freqs.items()),"maximum",max(freqs))
    print("Eligible offsets through 63:",[i for i in range(2,min(limit,63)+1) if candidate[i]])


if __name__=="__main__":
    parser=ArgumentParser(description=__doc__)
    parser.add_argument("--limit",type=int,default=2**20)
    args=parser.parse_args()
    print(dict(carry_checks()))
    print(dict(exact_pairs()))
    offset_counts(args.limit)
