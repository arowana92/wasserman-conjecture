"""Check the spread inequalities on exact coefficients and the prime-power family by digit allocation."""

from itertools import combinations, combinations_with_replacement, product
from math import gcd, prod

from check_endpoint_pair_bound import allocate, prime_powers, rows
from check_three_factors import valuation


def exact_pairs(limit):
    distinct = cells = forcing = 0
    for n in range(3,limit+1):
        table = list(rows(n))
        weight = tuple(prime_powers(n))
        pred = tuple(prime_powers(n-1))
        for row,first in table:
            d,f = row[1]-row[0],row[2]-row[1]
            for _,second in table:
                common = gcd(first,second)
                D,J = gcd(n,common),gcd(n-1,common)
                U,V,_ = allocate(weight,first,second)
                S,T,_ = allocate(pred,first,second)
                assert D*U*V == n and J*S*T == n-1
                bins = [gcd(S,gcd(row[(i+1)%3],row[(i+2)%3])) for i in range(3)]
                assert prod(bins)==S
                for i in range(3):
                    G = gcd(row[(i+1)%3],row[(i+2)%3])
                    assert row[i]*G % (S*U**2) == 0
                    assert 27*(n-1)**2*U <= (row[i]*G*J)**2*D**3
                    cells += 1
                if d==0 or f==0:
                    continue
                g = gcd(d,f)
                r,s = d//g,f//g
                h = r+s
                a,b,c = r//bins[2],s//bins[0],h//bins[1]
                P = a*b*c
                assert P*S == r*s*h
                B,br = divmod(b*(g*(r+h)-1),bins[1]*bins[2])
                C,cr = divmod(a*(g*(s+h)+1),bins[0]*bins[1])
                assert br == cr == 0 and min(B,C)>0
                assert B*r+C*s==3*P*g
                assert B*r*(r+2*s)-C*s*(2*r+s)==-3*P
                assert (B-C)%c == 0
                assert 27*(n-1)**4*P*(g//U)**6 <= D**6*J**4*(d+f)*(d*f)**4
                distinct += 1
                forcing += 27*(n-1)**4 > (d+f)*(d*f)**4
    print('Exact coefficient tests:',cells,'part and gcd bounds;',distinct,
          'distinct-row pairs;',forcing,'force a common divisor by spread alone')


def distributed_rows(prime,m,subset):
    Q = prime**(16*m)
    units = [1,prime**(12*m)] + [prime**i for i in subset]
    # Use all ordered allocations of lower digits and unordered allocations of
    # the three indistinguishable leading units to the three labelled positions.
    for upper in combinations_with_replacement(range(3),3):
        start = [upper.count(i)*Q for i in range(3)]
        for destinations in product(range(3),repeat=len(units)):
            row = start.copy()
            for unit,i in zip(units,destinations):
                row[i] += unit
            if min(row)>0:
                yield sorted(row)


def check_family():
    count = weights = maximum = 0
    for p in (5,7):
        for m in range(1,6):
            for length in range(m+1):
                for subset in combinations(range(1,m+1),length):
                    D = 1+sum(p**i for i in subset)
                    c = p**(12*m)+D
                    n = 3*p**(16*m)+c
                    assert D**4*c**5 <= 27*(n-1)**4
                    maximum = max(n,maximum)
                    weights += 1
                    for row in distributed_rows(p,m,subset):
                        assert sum(row) == n
                        count += 1
                        if gcd(*row)<3:
                            continue
                        if any(row[i]*gcd(row[(i+1)%3],row[(i+2)%3])<9*(n-1) for i in range(3)):
                            continue
                        d,f = row[1]-row[0],row[2]-row[1]
                        assert d*f>0 and (d+f)*(d*f)**4 < 27*(n-1)**4
    print('Digit allocations:',count,'positive labelled rows at',weights,'weights; largest',maximum)


def check_four_primes():
    from check_scaled_row_blocks import certify_prime

    p=11
    n=3*p**32+p**24+p**2+1
    q=3604276009509249793981575479
    s=72204423293878581162091272241
    r=n-109
    weight={2:1,3:2,89:1,1097:1,q:1}
    pred={5:2,11:2,29:1,s:1}
    assert prod(k**v for k,v in weight.items())==n
    assert prod(k**v for k,v in pred.items())==n-1
    for prime in set(weight)|set(pred)|{r}:
        certify_prime(prime)
    free=set(map(tuple,distributed_rows(p,2,(2,))))
    assert len(free)==30
    assert all(valuation(n,row,base)>0 for row in free for base in (q,s,r))
    assert q*(2*s+1)>n and n-r < min(2*q,s) and 2*r>n
    for base in set(weight)|set(pred):
        x=n
        digits=0
        while x:
            x,digit=divmod(x,base)
            digits+=digit
        assert digits>=3
    print('Optimal four-prime certificate at',n,':',p,q,s,r)


if __name__=='__main__':
    exact_pairs(40)
    check_family()
    check_four_primes()
