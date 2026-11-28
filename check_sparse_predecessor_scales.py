"""Direct factorial-floor checks for exact sparse predecessor losses."""

from math import comb, gcd

from check_direction_pair_count import valuation
from check_direction_predecessor_loss import vp,constant


def defect(row,v,p):
    n=sum(row)*v
    return max(0,vp(n-1,p)-valuation(n,tuple(x*v for x in row),p))


def loss_236(v):
    d=1
    binary=3*v-1
    if binary>=8 and binary&(binary-1)==0 and (binary.bit_length()-1)%2:
        d=2
    ternary=2*v-1
    if ternary>=9 and 3**vp(ternary,3)==ternary:
        assert d==1
        d=3
    return d


def run():
    row=(2,3,6)
    vs=set(range(1,100001))
    vs.update((2**k+1)//3 for k in range(3,1600,2))
    vs.update((3**k+1)//2 for k in range(2,1200))
    for v in vs:
        h=2**defect(row,v,2)*3**defect(row,v,3)
        assert h==loss_236(v)
        if v<=40:
            n=11*v
            exact=comb(n,2*v)*comb(9*v,3*v)
            assert (n-1)//gcd(n-1,exact)==h
    print("2:3:6 exact loss checks",len(vs),flush=True)
    count=0
    for p in (2,3,5,7):
        for h in range(1,4):
            for m in range(h,4*h+1,h):
                b0=(p**m-1)//(p**h-1)
                for k in range(1,6):
                    b=p**k*b0
                    y=(1,b,(p**h-1)*b)
                    assert defect(y,1,p)==max(0,k+h-m)
                    for v in range(1,51):
                        assert defect(y,v,p)<=max(0,k+h-m)
                        count+=1
    # Directly compare the ray algorithm on a smaller subset.
    for p in (2,3,5):
        for m in range(1,5):
            for k in range(1,4):
                b=p**k*((p**m-1)//(p-1))
                y=(1,b,(p-1)*b)
                assert vp(constant(y),p)==max(0,k+1-m)
    print("prime-power digit-sum bounds",count)
    for k in range(1,7):
        y=(1,2**k,2**k)
        for v in range(1,4096):
            j=defect(y,v,2)
            if j:
                assert v%2==1 and bin(v).count('1')<=k-j+1
            assert (j==k)==(v==1)


if __name__=="__main__":
    run()
