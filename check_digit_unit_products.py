"""Digit extrema for a single row, a global product test and an optimal prime list."""

from argparse import ArgumentParser
from itertools import combinations
from math import comb, gcd, prod

from flint import fmpz

from check_scaled_row_blocks import certify_prime
from search_endpoint import smallest_factors


def distribute(loads, digit):
    """Add digit indistinguishable units to least-loaded positions in batches."""
    low, high = min(loads), min(loads)+digit+1
    while low+1 < high:
        middle = (low+high)//2
        if sum(max(0,middle-x) for x in loads)<=digit:
            low = middle
        else:
            high = middle
    new = [max(low,x) for x in loads]
    left = digit-sum(b-a for a,b in zip(loads,new))
    for i,x in enumerate(new):
        if left and x==low:
            new[i] += 1
            left -= 1
    assert left == 0
    return new


def extrema(n, p):
    digits = []
    lowest = []
    weight, power = n, 1
    while weight:
        weight, digit = divmod(weight,p)
        digits.append(digit)
        lowest.extend([power]*min(digit,3-len(lowest)))
        power *= p
    if len(lowest)<3:
        return None
    loads = [0,0,0]
    for digit in reversed(digits):
        loads = distribute([p*x for x in loads],digit)
    assert min(loads)>0 and sum(loads)==n
    a,b,_ = lowest
    return prod(lowest), a*b*(n-a-b), prod(loads)


def exact_checks():
    top = 120
    spf = smallest_factors(top)
    tests = occupied = universal = 0
    for n in range(3,top+1):
        rows = [(a*b*(n-a-b),comb(n,a)*comb(n-a,b)) for a in range(1,n//3+1)
                for b in range(a,(n-a)//2+1)]
        for p in range(2,n+1):
            if spf[p]!=p:
                continue
            tests += len(rows)
            values = [m for m,k in rows if k%p]
            computed = extrema(n,p)
            if not values:
                universal += 1
                assert computed is None
                continue
            occupied += 1
            g = gcd(*values)
            power = 1
            while g%p==0:
                g//=p
                power*=p
            assert computed==(power,min(values),max(values)), (n,p,computed)
    print('exact coefficient tests',tests,'nonempty extrema',occupied,'universal cases',universal)


def four_primes():
    p,q = 1000000000039,1000000000193
    n = p*q
    bases = p,q,(n-1)//2,n-400
    bounds = {}
    for r in bases:
        assert fmpz(r).is_prime()
        certify_prime(r)
        bounds[r] = extrema(n,r)
        assert bounds[r] is not None
    for r,s in combinations(bases,2):
        div = bounds[r][0]*bounds[s][0]
        lo,hi = max(bounds[r][1],bounds[s][1]),min(bounds[r][2],bounds[s][2])
        assert ((lo+div-1)//div)*div>hi
        assert div>hi or lo>hi
    assert n==2*bases[2]+1 and p!=q and min(p,q)>2
    print('optimal certificate',n,bases,'pair exclusions',6)


def low_factor(n, p):
    full = 1
    h = 0 if n%p==0 else 1
    assert (n-h)%p==0
    cofactor = n-h
    while cofactor%p==0:
        full*=p
        cofactor//=p
    ex = extrema(n,p)
    if ex is None:
        return None
    baseline = full**(3-h)
    assert ex[0]%baseline==0
    return ex[0]//baseline


def controls():
    counts = [0,0,0]
    for p in (2,3,5,7,11):
        for a in range(1,9):
            for b in range(a,a+6):
                for c in (b+1,b+4,b+17):
                    for t in (1,7,799*p**(b-a),800*p**(b-a)):
                        n=p**a*(1+p**b+t*p**c)
                        w=low_factor(n,p)
                        assert w is None or w>=p**(b+c)
                        assert n+3<=864*p**(b+c)
                        counts[0]+=1
    # CRT classes with one low unit in both reduced bases.
    for b in range(1,21):
        u,v=2**(b+1),3**(b+1)
        residue=2+u*(pow(u,-1,v)%v)
        last=(864*6**(2*b)-3-residue)//(u*v)
        for j in {0,1,last//2,last}:
            n=residue+u*v*j
            assert n%u==2 and n%v==3 and n+3<=864*6**(2*b)
            x,y=low_factor(n,2),low_factor(n,3)
            assert x is None or y is None or x*y>=6**(2*b)
            counts[1]+=1
    # Every possible positive row maximum in base two has the explicit closed form.
    for n in range(3,3000):
        bits=[1<<k for k in range(n.bit_length()) if n & (1<<k)]
        ex=extrema(n,2)
        if len(bits)>=3:
            largest,second=bits[-1],bits[-2]
            assert ex[2]==largest*second*(n-largest-second)
            counts[2]+=1
    print('three-block, simultaneous congruence and binary controls',*counts)


def admits_partition(bounds, maximum):
    def extend(state, entry):
        divisor,low,high = state
        d,l,h = entry
        divisor *= d
        low,high = max(low,l),min(high,h)
        if (low+divisor-1)//divisor*divisor<=high:
            return divisor,low,high
        return None

    def visit(i, left, right):
        if i==len(bounds):
            return True
        nxt = extend(left,bounds[i])
        if nxt is not None and visit(i+1,nxt,right):
            return True
        if left==right:
            return False
        nxt = extend(right,bounds[i])
        return nxt is not None and visit(i+1,left,nxt)

    empty = (1,1,maximum)
    return visit(0,empty,empty)


def survey():
    counts = [0]*4
    for n in range(10**12+1,10**12+10001):
        factors = {int(p) for v in (n-1,n) for p,_ in fmpz(v).factor()}
        entries = {p:extrema(n,p) for p in factors}
        if None in entries.values():
            counts[0]+=1
            continue
        if 27*prod(x[0] for x in entries.values())>n**3*entries[2][2]:
            counts[1]+=1
            continue
        # Most restrictive required divisors first.
        if not admits_partition(sorted(entries.values(),reverse=True),n**3//27):
            counts[2]+=1
        else:
            counts[3]+=1
    print('ten-thousand-weight survey (universal, product, partition, unresolved)',*counts)


def survey_crt():
    from check_power_partner_assignments import check_weight
    weights = assignments = rows = 0
    for n in range(10**12+1,10**12+10001):
        count, found = check_weight(n)
        weights += count>0
        assignments += count
        rows += found
        if found:
            raise ArithmeticError(('inspect surviving row assignment',n))
    print('independent CRT audit, weights with assignments / assignments / rows',weights,assignments,rows)


if __name__=='__main__':
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('--crt',action='store_true')
    args=parser.parse_args()
    exact_checks()
    four_primes()
    controls()
    survey()
    if args.crt:
        survey_crt()
