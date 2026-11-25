"""Check polynomial specialization against exact coefficients and interval coverage by sparse parts."""

from collections import Counter
from itertools import combinations, product
from math import comb, gcd, prod

from sympy import Poly, symbols

from check_polynomial_powers import assignment_bounds, divisor_bound, evaluate, full_power_bound, resultant
from search_endpoint import atoms, smallest_factors


def compare_resultants():
    x = symbols("x")
    polys = [(a,b,c) for a,b,c in product(range(-1,3),range(3),range(1,3))]
    polys += [(1,1),(0,1),(2,3),(3,)]
    exact = {f:Poly.from_list(list(reversed(f)),gens=x) for f in polys}
    for f,g in product(polys,repeat=2):
        if resultant(f,g) != abs(exact[f].resultant(exact[g])):
            raise ArithmeticError(("resultant",f,g))
    # Each pair shares a polynomial factor, but all three together are coprime.
    inputs = ((0,1,1),(0,2,1),(2,3,1))
    b, = divisor_bound(*inputs)
    if not b or any(b % gcd(*(evaluate(f,t) for f in inputs)) for t in range(1,100)):
        raise ArithmeticError(("multi-polynomial divisor",b))
    # A proper power in the bound is missing part of the full power in the target.
    F,L = (2,5,4,1),(6,9,3)  # (x+1)^2(x+2), 3(x+1)(x+2)
    refined = full_power_bound(F,L)
    if len(refined) != 2:
        raise ArithmeticError(("unitary degree", refined))
    spf = smallest_factors(evaluate(F,60))
    for t in range(1,61):
        full = prod(w for _,w in atoms(evaluate(F,t),spf) if evaluate(L,t) % w == 0)
        if evaluate(refined,t) % full:
            raise ArithmeticError(("lost a unitary divisor",t,full))
    print(f"{len(polys)**2} independent resultants and 60 unitary specialization checks agree.")


def exact_rows():
    counts = Counter()
    spf = smallest_factors(900)
    for G in ((1,1,1),(2,1,1),(2,3,1),(1,0,2,1),(1,1,1,1),(2,2,2,1),(1,2,1,1,1)):
        for P,p in ((2,2),(3,3),(4,2),(5,5),(7,7)):
            Q = evaluate(G,P)
            if P <= max(G) or P*Q+1 >= len(spf):
                continue
            for eps in (0,1):
                n = P*Q+eps
                for a in range(1, Q//2+1):
                    for b in range(a,Q-a):
                        c = Q-a-b
                        if not eps and b > c:
                            continue
                        row = (P*a,P*b,P*c+eps)
                        coefficient = comb(n,row[0])*comb(n-row[0],row[1])
                        if coefficient % p == 0:
                            continue
                        forms = []
                        for part in (a,b,c):
                            digits = []
                            while part:
                                part,digit = divmod(part,P)
                                digits.append(digit)
                            forms.append(tuple(digits))
                        first,second = assignment_bounds(G,tuple(forms),2 if eps else None)
                        u = prod(w for r,w in atoms(n,spf) if coefficient % r)
                        s = prod(w for r,w in atoms(n-1,spf) if coefficient % r)
                        u,s = (u,s//P) if eps else (u//P,s)
                        if evaluate(first,P) % u or evaluate(second,P) % s:
                            raise ArithmeticError(("exact coefficient",n,row,u,s,first,second))
                        counts[eps] += 1
                        counts["nonconstant bound"] += max(len(first),len(second)) > 1
    if not counts["nonconstant bound"]:
        raise ArithmeticError("missing nonconstant control")
    print(f"Maximal assignments agree for exact coefficients: {dict(counts)}.")


def independent_interval():
    counts = Counter()
    for exponent in range(1,63):
        P = 2**exponent
        base = P*(P*(P*(P+1024)+1)+1)
        for offset in range(2,1021):
            n = base+offset
            if n <= 10**9:
                continue
            counts["weights"] += 1
            digits = [1 << i for i in range(n.bit_length()) if n & (1 << i)]
            # At least one of three disjoint positive parts uses at most a third of the ones.
            for k in range(1,len(digits)//3+1):
                for subset in combinations(digits,k):
                    z = sum(subset)
                    if z <= 2:
                        continue
                    u = gcd(n-z,z)
                    # Consecutive factors are coprime, so separate their gcds.
                    s = gcd(n-1,z-1)*gcd(n-1,z)
                    counts["sparse parts"] += 1
                    if 27*(n-1)**2 <= u**3*s**2:
                        raise ArithmeticError(("eligible sparse part",n,z))
    if counts["weights"] != 57064:
        raise ArithmeticError(("interval coverage",counts))
    print(f"Independent sparse-part coverage: {dict(counts)}.")


if __name__ == "__main__":
    compare_resultants()
    exact_rows()
    independent_interval()
