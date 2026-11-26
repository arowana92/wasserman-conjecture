"""Independently recover the small-product assignments in increasing bin order."""

from argparse import ArgumentParser
from collections import Counter
from functools import lru_cache
from math import gcd, isqrt, prod
from multiprocessing import Pool

from check_extremal_parts import divisors, factors


@lru_cache(maxsize=32768)
def positive_divisors(n):
    return divisors(factors(n))


def layouts(n):
    # List an unordered triple of multipliers from its product.
    for a in range(1, n+1):
        if a**3 > n:
            break
        if n % a:
            continue
        for b in range(a, isqrt(n//a)+1):
            if n % (a*b) == 0:
                yield (a,b,n//(a*b))


def build(left, right, direct=False, shard=0, shards=1):
    result = set()
    full = (left, right)

    def visit(rows, pending, bound):
        count = sum(map(len,pending))
        if count == 3 and shards > 1 and hash(tuple(rows)) % shards != shard:
            return
        denominators = [prod(b for a,b in row) for row in rows]
        numer = [sum(a*(mod//b) for a,b in row) for row,mod in zip(rows,denominators)]
        S,T = denominators
        X,Y = numer
        D = S*T-X*Y
        if D <= 0:
            return
        if count == 1:
            side = int(not pending[0])
            coefficient = pending[side][0]*denominators[side]*numer[1-side]
            v, rem = divmod(coefficient-1,D)
            if rem or v < bound or gcd(v, S*T) != 1:
                return
            final = [list(r) for r in rows]
            final[side].append((pending[side][0],v))
            moduli = [prod(b for a,b in row) for row in final]
            scales = [sum(a*(mod//b) for a,b in row) for row,mod in zip(final,moduli)]
            assert prod(scales) == prod(moduli)+1
            if gcd(*scales) == 1:
                result.add(tuple(sorted(tuple(sorted(row)) for row in final)))
            return

        # The next bin v is the smallest still missing, on either side.
        # Upper-bound each unknown fraction by its numerator divided by v.
        A,B = map(sum,pending)
        linear = X*T*B + Y*S*A
        cap = (linear+isqrt(linear*linear+4*D*S*T*A*B))//(2*D)
        for side in (0,1):
            for c in set(pending[side]):
                remaining = list(pending)
                unused = list(remaining[side])
                unused.remove(c)
                remaining[side] = tuple(unused)
                candidates = range(bound, cap+1)
                if count == 2 and cap-bound > 96 and not direct:
                    # Fit the bilinear equation by evaluating its cleared polynomial.
                    last_side = int(not remaining[0])
                    last = remaining[last_side][0]

                    def polynomial(v, w):
                        complete = list(rows)
                        complete[side] += ((c,v),)
                        complete[last_side] += ((last,w),)
                        sides = []
                        bins = []
                        for entries in complete:
                            (f,x),(g,y),(h,z) = entries
                            sides.append(f*y*z+g*x*z+h*x*y)
                            bins.append(x*y*z)
                        return prod(bins)-prod(sides)+1

                    constant = polynomial(0,0)
                    linear_v = polynomial(1,0)-constant
                    linear_w = polynomial(0,1)-constant
                    kernel = linear_v*linear_w-D*constant
                    if kernel <= 0:
                        continue
                    candidates = [(r-linear_w)//D for r in positive_divisors(kernel)
                                  if (r-linear_w)%D == 0 and bound*D <= r-linear_w <= cap*D]
                for v in candidates:
                    if gcd(v, S*T) == 1:
                        new_rows = list(rows)
                        new_rows[side] += ((c,v),)
                        visit(new_rows, remaining, v)

    visit([(),()], list(full), 1)
    return result


def digit_sum(n, prime):
    total = 0
    while n:
        n, digit = divmod(n, prime)
        total += digit
    return total


def audit(n, direct, jobs=1):
    tasks = []
    for d in range(1, isqrt(n)+1):
        if n % d == 0:
            for f in layouts(d):
                for g in layouts(n//d):
                    tasks.extend((f,g,direct,shard,jobs) for shard in range(jobs))
    found = set()
    if jobs == 1:
        for task in tasks:
            found.update(build(*task))
    else:
        with Pool(jobs) as pool:
            for systems in pool.starmap(build,tasks,chunksize=1):
                found.update(systems)
    primes = Counter()
    largest = 0
    for system in found:
        values = [[c*prod(v for k,(_,v) in enumerate(row) if k != j)
                   for j,(c,_) in enumerate(row)] for row in system]
        weight = prod(map(sum,values))
        largest = max(largest,weight)
        rows = [[z*sum(values[1-side]) for z in row] for side,row in enumerate(values)]
        if all(any(a&b for i,a in enumerate(row) for b in row[i+1:]) for row in rows):
            continue
        for q,_ in sorted(factors(weight)+factors(weight-1)):
            if all(sum(digit_sum(z,q) for z in row)>digit_sum(weight,q) for row in rows):
                primes[q] += 1
                break
        else:
            raise ArithmeticError(('check candidate',weight,rows))
    print(n, len(found), 'largest weight', largest, 'common primes after binary', dict(primes),
          flush=True)


def construction(limit):
    for n in range(2,limit+1,42):
        a=n+1
        b=n*a+1
        c=n*a*b+1
        v=a*b+a*c+b*c
        d=v+1
        e=v*d-1
        u=n*d*e+d+e
        weight=u*v
        bins=(a,b,c,d,e)
        assert all(gcd(x,y)==1 for i,x in enumerate(bins) for y in bins[i+1:])
        assert gcd(u,v)==1 and weight-prod(bins)==1
        x=[u*b*c,u*a*c,u*a*b]
        y=[v*n*d*e,v*e,v*d]
        assert sum(x)==sum(y)==weight
        assert prod(x)*prod(y)==n*weight**3*(weight-1)**2
        assert (e+d)%v==0 and gcd(e,d)==1 and gcd(*x)==u
    print(f'Checked the construction at every parameter 2 modulo 42 through {limit}.')


if __name__ == '__main__':
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('limit',type=int,nargs='?',default=32)
    parser.add_argument('--direct',action='store_true',help='scan every fifth bin directly')
    parser.add_argument('--jobs',type=int,default=1,help='number of processes for the prefix search')
    parser.add_argument('--start',type=int,default=2,help='first even product to check')
    args=parser.parse_args()
    if args.jobs < 1 or args.start < 2 or args.start % 2:
        parser.error('use positive jobs and an even positive starting product')
    construction(42_002)
    for n in range(args.start,args.limit+1,2):
        audit(n,args.direct,args.jobs)
