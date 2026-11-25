"""Certify trinomial weights P*G(P) and P*G(P)+1 for bounded polynomial coefficients."""

from argparse import ArgumentParser
from collections import Counter
from functools import lru_cache
from itertools import combinations, product
from math import gcd, lcm

from flint import fmpz_mat, fmpz_poly

from check_close_powers import INITIAL_LIMIT, bases, check_complement, prime_base, subproducts
from search_even_endpoint import avoids, prime_powers


def trim(coefficients):
    values = list(coefficients)
    while values and values[-1] == 0:
        values.pop()
    return tuple(values)


def evaluate(poly, x):
    value = 0
    for c in reversed(poly):
        value = value*x+c
    return value


def row_shapes(G):
    # The monic leading unit fixes slot 3. The other two slots are unordered.
    local = [[(i,j,k-i-j) for i in range(k+1) for j in range(k-i+1)] for k in G[:-1]]
    for blocks in product(*local):
        a,b,c = (trim(k) for k in zip(*blocks,(0,0,1)))
        if a and b and a <= b:
            yield a,b,c


@lru_cache(maxsize=100000)
def resultant(f,g):
    """Absolute Sylvester determinant for nonzero coefficient tuples in ascending order."""
    m,n = len(f)-1,len(g)-1
    if min(m,n) == 0:
        return abs(f[0])**n if m == 0 else abs(g[0])**m
    rows = [[0]*i+list(f)+[0]*(n-1-i) for i in range(n)]
    rows += [[0]*i+list(g)+[0]*(m-1-i) for i in range(m)]
    return abs(int(fmpz_mat(rows).det()))


def cap(*polys):
    # Every common evaluation divisor divides each constant and each
    # nonzero pair resultant. Zero denotes the absence of a finite bound.
    bound = gcd(*(abs(f[0]) for f in polys if len(f) == 1))
    for f,g in combinations(polys,2):
        if bound == 1:
            break
        if len(f) > 1 and len(g) > 1 and f != g:
            bound = gcd(bound,resultant(f,g))
    return bound


def coeffs(poly):
    return tuple(map(int,poly.coeffs()))


@lru_cache(maxsize=65536)
def common_bound(polys):
    bound = cap(*polys)
    if bound:
        return (bound,)
    fs = [fmpz_poly(list(f)) for f in polys]
    shared = fs[0]
    for g in fs[1:]:
        shared = shared.gcd(g)
    reduced = [g//shared for g in fs]
    bound = cap(*(coeffs(g) for g in reduced))
    if not bound:
        # Some pairwise gcds can be nonconstant even when the whole list is coprime.
        # For each root of reduced[0] only finitely many integer t make the sum vanish.
        t = 1
        while not bound:
            mixed = sum((t**j*g for j,g in enumerate(reduced[1:])),fmpz_poly([]))
            if mixed:
                bound = cap(coeffs(reduced[0]),coeffs(mixed))
            t += 1
    return coeffs(bound*shared)


def divisor_bound(*polys):
    return common_bound(tuple(sorted(set(polys))))


@lru_cache(maxsize=32768)
def common_multiple(polys):
    value = fmpz_poly([1])
    for f in set(polys):
        g = fmpz_poly(list(f))
        value = (value//value.gcd(g))*g
    return coeffs(value)


def assignment_bounds(G, row, shifted):
    if shifted is None:
        m = (-1,)+G
        v = divisor_bound(G,*row)
        s = common_multiple(tuple(divisor_bound(m,f,g) for f,g in combinations(row,2)))
        return v,s
    n = (1,)+G
    u = divisor_bound(n,*(row[i] for i in range(3) if i != shifted))
    groups = []
    for k in range(3):
        terms = [G]
        for i,part in enumerate(row):
            delta = int(i == shifted)-int(i == k)
            terms.append((delta,)+part if delta else part)
        groups.append(divisor_bound(*terms))
    return u,common_multiple(tuple(groups))


@lru_cache(maxsize=65536)
def power_limit(G, eps, first, second):
    if eps and len(first) == 1 and first[0] < 3:
        return 0
    n = len(G)-1
    v = sum(map(abs,first))
    s = sum(map(abs,second))
    exponent = 2*n - (1-eps) - 3*(len(first)-1) - 2*(len(second)-1)
    # A quick strict degree bound often covers the whole initial interval.
    if exponent > 0:
        bound = (1 if eps else 4)*v**3*s**2//27
        hi = 1
        while hi**exponent <= bound:
            hi *= 2
        if hi*evaluate(G,hi)+eps <= INITIAL_LIMIT:
            return hi
    x = fmpz_poly([0,1])
    q,f,g = [fmpz_poly(list(poly)) for poly in (G,first,second)]
    # Square so signs at a particular evaluation are immaterial.
    diff = 729*q**4-f**6*g**4 if eps else 729*(x*q-1)**4-x**6*f**6*g**4
    if not diff or diff[diff.degree()] < 0:
        return None
    # At a shift with positive constant and nonnegative Taylor coefficients
    # the necessary inequality fails on the entire half-line.
    def positive(k):
        shifted = diff(x+k)
        return shifted[0] > 0 and all(c >= 0 for c in shifted.coeffs())
    if positive(0):
        return 0
    lo,hi = 0,1
    while not positive(hi):
        lo,hi = hi,2*hi
    while hi-lo > 1:
        mid = (lo+hi)//2
        if positive(mid):
            hi = mid
        else:
            lo = mid
    return hi-1


@lru_cache(maxsize=16384)
def full_power_bound(target, bound):
    # The assigned divisor is unitary in target(P). A partial polynomial
    # factor in the bound can therefore contribute only a bounded integer.
    if len(bound) == 1:
        return bound
    F = fmpz_poly(list(target))
    L = fmpz_poly(list(bound))
    content, factors = L.factor()
    result = fmpz_poly([content])
    for f,k in factors:
        quotient = F
        for _ in range(k):
            quotient, remainder = divmod(quotient,f)
            if remainder:
                raise ArithmeticError(("bound factor outside target",target,bound))
        if quotient % f:
            result *= f**k
        else:
            rest = L//(f**k)
            c = cap(coeffs(f),coeffs(rest))
            if not c:
                raise ArithmeticError(("noncoprime remainder",f,rest))
            result *= c**k
    return coeffs(result)


def test_assignment(G, row, shifted, P, counts):
    Q = evaluate(G,P)
    n = P*Q + int(shifted is not None)
    entries = tuple(P*evaluate(v,P)+int(i == shifted) for i,v in enumerate(row))
    u = gcd(*entries)
    s = gcd(n-1,lcm(*(gcd(f,g) for f,g in combinations(entries,2))))
    if u < 3 or u**3*s*s < 27*(n-1)**2:
        return
    counts["above numerical size bound"] += 1
    p = prime_base(P)
    if not p or not avoids(n,entries,p):
        return
    # Take all eligible full powers for this first row. If its complement
    # worked for a smaller assignment, it works after removing obligations.
    U,S = (1,1)
    for r,w in prime_powers(n).items():
        if all(v % w == 0 for v in entries) and avoids(n,entries,r):
            U *= w
    if U < 3 or n//U < 3:
        return
    for r,w in prime_powers(n-1).items():
        if sorted(v % w for v in entries) == [0,0,1] and avoids(n,entries,r):
            S *= w
    T = (n-1)//S
    if 27*T*T <= U**3:
        counts["complement tests"] += 1
        check_complement(n,U,T,counts)


def irreducible(f):
    _, factors = fmpz_poly(f).factor()
    return len(factors) == 1 and factors[0][1] == 1


def check_box(degree,height,list_exceptions=False):
    counts = Counter()
    for coeff in product(range(height+1),repeat=degree):
        G = coeff+(1,)
        if height*evaluate(G,height)+1 > INITIAL_LIMIT:
            raise ValueError("P through the coefficient height must fit the verified initial interval")
        complete = [True,True]
        base_irreducible = irreducible(list(G))
        generic = [base_irreducible and irreducible([s,*G]) for s in (-1,1)]
        low = max(2,max(G)+1)
        for row in row_shapes(G):
            for shifted in (None,0,1,2):
                if shifted == 1 and row[0] == row[1]:
                    continue
                eps = int(shifted is not None)
                counts["shapes",eps] += 1
                first,second = assignment_bounds(G,row,shifted)
                high = power_limit(G,eps,first,second)
                if high is None:
                    targets = ((1,)+G,G) if eps else (G,(-1,)+G)
                    refined = tuple(full_power_bound(f,b) for f,b in zip(targets,(first,second)))
                    if refined != (first,second):
                        high = power_limit(G,eps,*refined)
                        if high is not None:
                            counts["resolved by unitarity",eps] += 1
                if high is None:
                    complete[eps] = False
                    if generic[eps]:
                        raise ArithmeticError(("unbounded irreducible family",G,row,shifted))
                    counts["unbounded shapes",eps] += 1
                    continue
                if high < low or high*evaluate(G,high)+eps <= INITIAL_LIMIT:
                    continue
                counts["above initial interval",eps] += 1
                for P in range(low,high+1):
                    if P*evaluate(G,P)+eps > INITIAL_LIMIT:
                        counts["integer evaluations",eps] += 1
                        test_assignment(G,row,shifted,P,counts)
        for eps in (0,1):
            counts["irreducible pairs",eps] += generic[eps]
            counts["complete",eps] += complete[eps]
        if list_exceptions and not all(complete):
            print(f"G coefficients {G}, unresolved offsets {[i for i in (0,1) if not complete[i]]}")
        counts["complete pairs"] += all(complete)
    print(f"Degree {degree}, height {height}: {dict(counts)}",flush=True)
    return counts


if __name__ == "__main__":
    cli = ArgumentParser(description=__doc__)
    cli.add_argument("degree",type=int)
    cli.add_argument("height",type=int)
    cli.add_argument("--list-exceptions",action="store_true")
    args = cli.parse_args()
    if min(args.degree,args.height) < 1:
        cli.error("use positive degree and height")
    check_box(args.degree,args.height,args.list_exceptions)
