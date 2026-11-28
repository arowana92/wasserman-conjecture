"""Explicit minimum valuations on (p-1)X+pY-p(p^k-1)Z = e."""

import argparse
from itertools import product
from math import comb

from check_line_divisors import minimum, prime_factors, valuation


def digit_sum(x, p):
    total = 0
    while x:
        x, digit = divmod(x,p)
        total += digit
    return total


def lowest_tokens(n, count, p):
    block = 0
    power = 1
    while count and n:
        n, digit = divmod(n,p)
        take = min(count,digit)
        block += take*power
        count -= take
        power *= p
    if count:
        raise ValueError("source has too few digits")
    return block


def high_power(n,p):
    power = 1
    while power <= n:
        power *= p
    return power


def minimize(p, k, offset):
    """Preconditions: p prime, k positive. Return exponent and a positive row."""
    if k<1 or p<2 or prime_factors(p)!=[p]:
        raise ValueError("p must be prime and k must be positive")
    u = -offset % p
    v = (offset+u)//p
    pk = p**k
    if v<0 or digit_sum(v,p)>u:
        if v<0:
            # Keep spare full digits below the digit occupied by z.
            z = p**3*high_power(-v,p)
            source = pk*z+v
        else:
            z = 0
            source = v
        block = lowest_tokens(source,u,p)
        t = (block-u)//(p-1)
        y = source-block-z
        exponent = 0
    elif u==p-1 and v==0:
        # Y=p-1 is the smallest positive multiple of p-1. T<p^k.
        # A congruence gives the same construction uniformly, including k=1.
        b = (pk-1)//(p-1)
        t = (-2)%b or b
        z = (2+t)//b
        y = p-1
        exponent = 2
    else:
        s = digit_sum(v,p)
        missing = u-s
        z = high_power(v,p)
        high = (pk//p)*z
        t = (v-s + missing*(high-1) + high-z)//(p-1)
        y = (p-1-missing)*high
        exponent = 1
    x = p*t+u
    if x==0 or z==0:
        shift = high_power(x+y+z,p)
        x += p*(pk-1)*shift
        z += (p-1)*shift
    row = (x,y,z)
    assert min(row)>0 and (p-1)*x+p*y-p*(pk-1)*z==offset
    assert valuation(row,p)==exponent
    return exponent, row


def odd_witness(k, offset, q):
    """A q-free positive row on the binary line, q an odd prime."""
    t = 2*(2**k-1)
    high = high_power(abs(offset),q)
    lead = t-int(offset<0)
    use_y = int(lead%q==q-1)
    row = ((t-2*use_y)*high+offset,use_y*high,high)
    if row[0]==0:
        assert (k,offset,q)==(1,0,3)
        row = (4,1,3)
    if row[1]==0:
        cap = q
        while (cap-1)%(2**k):
            cap *= q
        z = (cap-1)//(2**k)
        shift = high_power(sum(row),q)
        row = (row[0],shift*(2**k-1)*z,row[2]+shift*z)
    assert min(row)>0 and row[0]+2*row[1]-t*row[2]==offset
    assert valuation(row,q)==0
    return row


def audit():
    searches = exact = 0
    for p, max_k, radius in [(2,4,120),(3,3,90),(5,2,80),(7,1,80)]:
        for k, offset in product(range(1,max_k+1),range(-radius,radius+1)):
            exponent, row = minimize(p,k,offset)
            other, witness = minimum(p-1,p,p*(p**k-1),p,offset)
            assert exponent==other, (p,k,offset,row,witness)
            if sum(row)<=1500:
                n = sum(row)
                coefficient = comb(n,row[0])*comb(n-row[0],row[1])
                assert coefficient%p**exponent==0 and coefficient%p**(exponent+1)!=0
                exact += 1
            searches += 1
    large = 0
    for p, k, j in product((2,3,5,7,11), (1,2,5,8), (5,20,80)):
        for v in (-p**j-1,-p**j,-p**j+1,0,p**j-1,p**j,p**j+1):
            for u in range(p):
                minimize(p,k,p*v-u)
                large += 1
    odd = 0
    for k, q, offset in product(range(1,7),(3,5,7,11,13,17,19,23,29,31),range(-60,61)):
        odd_witness(k,offset,q)
        odd += 1
    print(f"odd-prime witnesses={odd}")
    print(f"shortest-path comparisons={searches}, direct factorial-ratio checks={exact}, large-offset checks={large}")


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--line",type=int,nargs=3,metavar=("P","K","OFFSET"))
    args = parser.parse_args()
    if args.line:
        print(minimize(*args.line))
    else:
        audit()
