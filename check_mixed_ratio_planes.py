"""All-weight mixed relation certificates using exact divisors and affine carries."""

import argparse
from collections import Counter, deque
from functools import lru_cache
from itertools import combinations, combinations_with_replacement, product
from math import gcd, lcm, prod

from check_affine_endpoint import divisors, factors
from carry_blocks import CarryGraph
from check_part_ratio_divisors import sharp_constant


@lru_cache(maxsize=None)
def pair_data(i, j):
    h = gcd(i, j)
    a, b = sorted((i//h, j//h))
    return sharp_constant(a, b), a+b


def passes_pair_bounds(cap, pair_sum, row):
    for i, j in combinations(row, 2):
        other_cap, other_sum = pair_data(i, j)
        if lcm(cap, other_cap) <= pair_sum*other_sum:
            return False
    return True


@lru_cache(maxsize=None)
def divisors_of(n):
    return tuple(divisors(factors(n)))


def allocations(cap, normal):
    """Assign each full power of D to the ratio or a coefficient of the plane."""
    for groups in product(*(divisors_of(abs(v)) for v in (cap, *normal))):
        d = prod(groups)
        if d == lcm(*groups):
            yield d, groups[1:]


def plane_rows(normal, steps, limit):
    """Choose two multiples, then solve the homogeneous equation for the third."""
    i, j, k = sorted(range(3), key=lambda pos: -steps[pos])
    row = [0]*3
    for x in range(steps[i], limit-steps[j]-steps[k]+1, steps[i]):
        for y in range(steps[j], limit-x-steps[k]+1, steps[j]):
            z, rem = divmod(-normal[i]*x-normal[j]*y, normal[k])
            if rem or z < 1 or z % steps[k] or x+y+z > limit:
                continue
            row[i], row[j], row[k] = x, y, z
            yield tuple(row)


def shapes(height):
    found = set()
    positive = range(1, height+1)
    ratios = [(a,b) for a,b in combinations_with_replacement(positive,2) if gcd(a,b)==1]
    for r,s in combinations_with_replacement(positive, 2):
        for t in positive:
            if gcd(r,s,t) != 1:
                continue
            normal = (r,s,-t)
            for a,b in ratios:
                cap = sharp_constant(a,b)
                for d, groups in allocations(cap, normal):
                    total = prod(groups)
                    steps = tuple(total//v for v in groups)
                    limit = (d-1)//(a+b)
                    if sum(steps) > limit:
                        continue
                    for row in plane_rows(normal, steps, limit):
                        u = sum(row)
                        if gcd(u,d)==1 and passes_pair_bounds(cap,a+b,row):
                            found.add((a,b,d,u,*sorted(row)))
    return found


def group_shapes(height):
    """Allow each of the three assigned companion products to reach height."""
    found = set()
    positive = range(1,height+1)
    for a,b in combinations_with_replacement(positive,2):
        if gcd(a,b)>1:
            continue
        cap = sharp_constant(a,b)
        for groups in combinations_with_replacement(positive,3):
            total = prod(groups)
            if total!=lcm(*groups):
                continue
            k,l,m = (total//v for v in groups)
            for r in divisors_of(cap):
                if gcd(total,r)>1:
                    continue
                d = r*total
                limit = (d-1)//(a+b)
                for x in range(k,limit-l-m+1,k):
                    for y in range(l,limit-x-m+1,l):
                        for z in range(m,limit-x-y+1,m):
                            u = x+y+z
                            row = (x,y,z)
                            if gcd(d,u)==1 and passes_pair_bounds(cap,a+b,row):
                                found.add((a,b,d,u,*sorted(row)))
    return found


def affine_rows(a,b,d,u,j,row):
    m = -pow(d,-1,u) % u
    v = (d*m+1)//u
    slopes = (a*u*j, b*u*j, d-(a+b)*u*j)
    first = (tuple(s*u for s in slopes), (slopes[0]*m, slopes[1]*m, slopes[2]*m+1))
    second = (tuple(i*d for i in row), tuple(i*v for i in row))
    return first, second, v


def admissible(slopes, offsets, p, allowed_first):
    """Existence of a nonnegative integer parameter giving a carry-free row."""
    def successor(state, digit):
        terms = tuple(c+s*digit for c,s in zip(state,slopes))
        if sum(t%p for t in terms)<p:
            return tuple(t//p for t in terms)
        return None

    initial = {successor(offsets, digit) for digit in allowed_first}
    initial.discard(None)
    queue = deque(initial)
    seen = set(initial)
    while queue:
        state = queue.popleft()
        if state == (0,0,0):
            return True
        for digit in range(p):
            dest = successor(state,digit)
            if dest is not None and dest not in seen:
                seen.add(dest)
                queue.append(dest)
    return False


def suffix_graph(p,slopes,offsets,w,r):
    """Apply all Legendre layers through w at once, then test reachable tails."""
    initial = tuple(k*r+t for k,t in zip(slopes,offsets))
    q = p
    while q<=w:
        if sum(initial)//q != sum(x//q for x in initial):
            return False
        q *= p
    graph = CarryGraph(p,slopes,tuple(x//w for x in initial))
    return bool(graph.accept)


def local_certificate(a,b,d,u,j,row,p):
    """Return true when every locally eligible row has a carry at this prime."""
    sl = (a*u*j,b*u*j,d-(a+b)*u*j)
    off = (0,0,1)
    zero = (0,0,0)
    if u%p==0:
        w = dict(factors(u))[p]
        m = -pow(d,-1,u)%u
        classes = {(m+u*k)%(w*p) for k in range(p) if ((d*(m+u*k)+1)//u)%p}
        return not any(suffix_graph(p,sl,off,w*p,r) for r in classes)
    if d%p==0:
        w = dict(factors(d))[p]
        if a%w==0 or b%w==0:
            h = b if a%w==0 else a
            if suffix_graph(p,sl,off,w,pow(h*u*j,-1,w)):
                return False
        if sum(v%w==0 for v in row)>=2:
            return not suffix_graph(p,row,zero,w,pow(u,-1,w))
        return True
    return not any(suffix_graph(p,coeff,bias,p,k)
                   for coeff,bias in ((sl,off),(sl,zero),(row,zero)) for k in range(1,p))


def classify(templates, prime_hint=None, full_parameter=False):
    excluded = Counter()
    survivors = []
    for a,b,d,u,*row in sorted(templates):
        for j in range(1,(d-1)//(u*(a+b))+1):
            if gcd(d,j)>1:
                continue
            A,B,v = affine_rows(a,b,d,u,j,row)
            d_powers = dict(factors(d))
            primes = sorted({2,*d_powers,*(p for p,w in factors(u))})
            if prime_hint is not None:
                hint = prime_hint(a,b,d,u,j,tuple(row),primes)
                if hint is None:
                    survivors.append((a,b,d,u,j,*row))
                    continue
                assert hint in primes
                primes.remove(hint)
                primes.insert(0,hint)
            for p in primes:
                if not full_parameter:
                    if local_certificate(a,b,d,u,j,tuple(row),p):
                        excluded[p] += 1
                        break
                    continue
                if u%p==0:
                    allowed = (i for i in range(p) if (v+d*i)%p)
                    possible = admissible(*A,p,allowed)
                else:
                    q = d_powers.get(p,1)
                    first = (a%q==0 or b%q==0)
                    second = sum(x%q==0 for x in row)>=2
                    m = -pow(d,-1,u) % u
                    digits = [i for i in range(p) if q==1 or (m+u*i)%p]
                    possible = ((first and admissible(*A,p,digits)) or
                                (second and admissible(*B,p,range(p))))
                if not possible:
                    excluded[p] += 1
                    break
            else:
                survivors.append((a,b,d,u,j,*row))
    return excluded, survivors


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--height",type=int,default=12)
    parser.add_argument("--groups",action="store_true",help="allow all bounded companion groups")
    parser.add_argument("--full-parameter",action="store_true",help="use the full common progression in the carry graph")
    args = parser.parse_args()
    candidates = group_shapes(args.height) if args.groups else shapes(args.height)
    print('starting shapes',len(candidates),flush=True)
    from check_ratio_groups_cross import prime_possible
    def hint(a,b,d,u,j,row,primes):
        return next((p for p in primes if not prime_possible(a,b,d,u,j,row,p)),None)
    excluded, survivors = classify(candidates,hint,args.full_parameter)
    print("height",args.height,"shapes",len(candidates),"excluded",dict(sorted(excluded.items())),
          "surviving families",len(survivors),flush=True)
    for item in survivors[:30]:
        print(item)
    if args.height==10:
        assert len(candidates)==(623 if args.groups else 531) and not survivors
    if args.height==12:
        exceptions = [(4,11,13860,673,1,63,70,540),(7,9,2772,167,1,11,44,112)]
        assert len(candidates)==(27387 if args.groups else 10332)
        assert survivors==(exceptions if args.groups else [])
