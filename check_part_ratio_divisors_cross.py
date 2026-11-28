"""Check the predicted exponent bounds by constructing rows attaining each one."""

from collections import Counter
from itertools import combinations_with_replacement
from math import comb, gcd

from check_part_ratio_divisors import minimum_decrement, sharp_constant
from check_primitive_row_scales import powers


def layers(parts,p):
    n=sum(parts)
    val=0
    while n:
        n//=p
        parts=[a//p for a in parts]
        val += n-sum(parts)
    return val


def transport_cap(a,b,p,k):
    if a%p==0:
        a,b=b,a
    q=p**k
    other=b//q
    num=a+b-other
    if num%a:
        return k
    scale=num//a
    h=0
    while scale%p==0:
        h+=1
        scale//=p
    return min(h,k) if scale==1 and h else k


def audit(height):
    count=Counter()
    for a,b in combinations_with_replacement(range(1,height+1),2):
        if gcd(a,b)!=1:
            continue
        constant=1
        for p,q in powers(a*b):
            k=0
            power=1
            while power!=q:
                power*=p
                k+=1
            exponent=transport_cap(a,b,p,k)
            constant*=p**exponent
            if exponent!=k:
                count['reduced prime exponents']+=1
            # Independent lower certificate: evaluate the whole positive
            # trinomial and compare the predecessor valuation at one weight.
            cost,g=minimum_decrement(a,b,p,q)
            n,level=1,0
            while n<=(a+b)*g:
                n*=p
                level+=1
            n+=1
            assert g%p and level-layers((a*g,b*g,n-(a+b)*g),p)==exponent
            count['prime witnesses']+=1
        assert sharp_constant(a,b)==constant
        count['ratios']+=1
        count['reduced constants']+=constant<a*b
        count['noninteger ratios with C at most sum']+=a>1 and constant<=a+b
    return count


def small():
    count=0
    for a,b in combinations_with_replacement(range(1,9),2):
        if gcd(a,b)>1:
            continue
        c=sharp_constant(a,b)
        for g in range(1,10):
            for z in range(1,30):
                n=(a+b)*g+z
                coef=comb(n,z)*comb((a+b)*g,a*g)
                assert c*g*coef%(n-1)==0
                count+=1
    return count


if __name__=="__main__":
    print("exact multinomials",small())
    print(dict(audit(512)))
