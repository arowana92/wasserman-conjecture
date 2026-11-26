"""Independently enumerate higher prime powers using the gcd of three exponent gaps."""

from argparse import ArgumentParser
from collections import Counter
from itertools import permutations
from math import gcd, isqrt, prod

from check_sparse_companions import factor, inspect
from check_three_factors import valuation
from search_endpoint import smallest_factors


def small_weights(p, P, cap):
    """Build the normalized cofactor from its digits, including a nonzero units digit."""
    place = 1
    while P*p*place <= cap:
        place *= p

    def visit(place, val, units_left, used):
        if place == 0:
            if len(used) >= 3 and val % p:
                yield P*val, used
            return
        for digit in range(min(p-1, units_left, (cap//P-val)//place)+1):
            yield from visit(place//p, val+digit*place, units_left-digit,
                             used + (P*place,)*digit)

    yield from visit(place, 0, 5, ())


def at_most_two(p, P, prior, systems):
    cap = 1+isqrt(4*P**3//27)
    if cap <= prior:
        return
    for n, units in small_weights(p,P,cap):
        if n <= prior or n//P < 3:
            continue
        if 27*(n-1)**2 <= P**3:
            systems.add((n,1))
        if n%4 != 3 or len(units) != 5:
            continue
        for a,b,c,d,e in set(permutations(units)):
            x,y = b+c,d+e
            if (a&x) == (a&y) == (x&y) == 0:
                systems.add((n,2))
                break


def larger_products(p, a, P, prior, systems, counts):
    terms = [p**i for i in range(3*a)]
    # The singleton exponent i is an even multiple of g, and the two
    # other gaps are odd multiples. At least one lower exponent equals a.
    for g in range(1,(3*a-1)//2+1):
        M = terms[g]+1
        cap = 1+isqrt(P**3*M*M//27)
        if cap <= prior:
            continue
        factors = factor(M)
        top = max(i for i,v in enumerate(terms) if v <= cap)
        for i in range(2*g,top+1,2*g):
            if i < a:
                continue
            for d in range(g,top-a+1,2*g):
                for e in range(d,top-a+1,2*g):
                    if gcd(i,d,e) != g:
                        continue
                    for j in range(a,top-d+1):
                        # When i and j both exceed a, the other lower exponent is fixed.
                        max_l = a if i > a and j > a else top-e
                        for l in range(a,max_l+1):
                            if d == e and j > l:
                                continue
                            row = (terms[i], terms[j]+terms[j+d], terms[l]+terms[l+e])
                            n = sum(row)
                            if n <= prior or n > cap or n//P % p == 0 or valuation(n,row,p):
                                continue
                            counts['candidate gap rows'] += 1
                            S = prod(w for q,w in factors if (n-1)%(q*w) != 0
                                     and valuation(n,row,q) == 0)
                            if S > 2 and 27*((n-1)//S)**2 <= P**3:
                                systems.add((n,S))


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--limit', type=int, default=200000000)
    parser.add_argument('--prior', type=int, default=10**9)
    args = parser.parse_args()
    assert 2 <= args.limit and 0 <= args.prior <= 10**9
    top = isqrt(args.limit)
    spf = smallest_factors(top)
    counts = Counter()
    for p in range(2,top+1):
        if spf[p] != p:
            continue
        a,P = 2,p*p
        while P <= args.limit:
            systems = set()
            at_most_two(p,P,args.prior,systems)
            larger_products(p,a,P,args.prior,systems,counts)
            for n,S in sorted(systems):
                T = (n-1)//S
                required = factor(T)
                if len(required) >= 3 and max(w for q,w in required)*2 < P:
                    print(f'CRT {p} {P} {n} {S}')
                    inspect(p,P,n,S,counts)
            a += 1
            P *= p
    print(dict(counts))


if __name__ == '__main__':
    main()
