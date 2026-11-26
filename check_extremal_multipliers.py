"""Enumerate all equality-shape triples at a given multiplier, independently of bin gaps."""

from argparse import ArgumentParser
from collections import Counter
from math import gcd

from check_extremal_parts import divisors, factors
from check_three_factors import valuation


def search(first, last, direct, emit):
    counts = Counter()
    largest = 0
    for ell in range(first+(first%2),last+1,2):
        for z in range(5,2*ell+1,4):
            d = ell+z
            if gcd(d,ell*(ell+1)) != 1:
                continue
            counts['factor pairs sought'] += 1
            b = ell*d
            kernel = b*b+z
            if direct:
                small = (-b) % z
                roots = range(small or z,b+1,z)
            else:
                roots = divisors(factors(kernel))
            for r in roots:
                if r>b or kernel % r or (r+b) % z:
                    continue
                e = (r+b)//z
                if e<=d or (e-ell)%4!=1 or gcd(e,ell+1)!=1:
                    continue
                f, rem = divmod(kernel//r+b,z)
                if rem or f<=e or (f-ell)%4!=1 or gcd(f,ell+1)!=1:
                    continue
                u = d*e+d*f+e*f
                t = d*e*f
                assert t==ell*u+1 and gcd(t,u-1)==1
                n = (u-1)*t+1
                v = n//u
                if gcd(v,u)!=1:
                    continue
                assert gcd(d,e)*gcd(d,f)*gcd(e,f) == 1 and u*v==n
                counts['unitary triples'] += 1
                largest = max(largest,n)
                binary_fit = u & (n-u) == 0
                assert binary_fit == (valuation(n,(u,n-u),2)==0)
                row = (v*d*e,v*d*f,v*e*f)
                assert sum(row)==n and valuation(n,row,2)>0
                if emit:
                    print('BINS',d,e,f,u,n)
                if binary_fit:
                    counts['odd binomial']+=1
                    print('ODD',ell,d,e,f,u,n)
    print(first,last,dict(counts),'largest',largest)


def main():
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('first',type=int)
    parser.add_argument('last',type=int)
    parser.add_argument('--direct',action='store_true')
    parser.add_argument('--emit',action='store_true')
    args=parser.parse_args()
    if args.first<2 or args.last<args.first:
        parser.error('supply a positive interval starting at two or above')
    search(args.first,args.last,args.direct,args.emit)


if __name__ == '__main__':
    main()
