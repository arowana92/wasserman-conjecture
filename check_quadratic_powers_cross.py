"""Test quadratic divisibility bounds on maximal assignments from exact coefficients."""

from collections import Counter
from math import comb, gcd, isqrt, prod

from check_quadratic_powers import factor_bounds
from search_endpoint import atoms, smallest_factors


def verify():
    counts = Counter()
    spf = smallest_factors(17*(7*17+6)+1)
    for P,p in ((2,2),(3,3),(4,2),(5,5),(7,7),(8,2),(9,3),(11,11),(13,13),(17,17)):
        for h in range(1,8):
            for d in range(1,min(P,7)):
                Q = h*P+d
                for eps in (0,1):
                    n = P*Q+eps
                    disc = d*d+(1-2*eps)*4*h
                    digit_one = any(v == 1 or atoms(v,spf) == [(p,v)] for v in (h,d))
                    if not digit_one and (gcd(h,d) > 2 or (disc >= 0 and isqrt(disc)**2 == disc)):
                        continue
                    # Directly enumerate the normalized parts; recover their
                    # high and low coefficients only after the exact p-test.
                    high_powers = atoms(n,spf)
                    other_powers = atoms(Q if eps else n-1,spf)
                    for x in range(1,Q//2+1):
                        for y in range(x,Q-x):
                            z = Q-x-y
                            if not eps and z < y:
                                continue
                            row = (P*x,P*y,P*z+eps)
                            K = comb(n,row[0])*comb(n-row[0],row[1])
                            if K % p == 0:
                                continue
                            form = tuple(divmod(v,P) for v in (x,y,z))
                            first,second = factor_bounds(h,d,form,eps)
                            U = prod(w for r,w in high_powers if K % r)
                            S = prod(w for r,w in other_powers if K % r)
                            assigned = U if eps else U//P
                            if first % assigned or second % S:
                                raise ArithmeticError(("divisibility",n,row,assigned,S,first,second))
                            counts["exact rows",eps] += 1
                            counts["nontrivial first bound",eps] += assigned > 1
                            counts["nontrivial second bound",eps] += S > 1
    if not all(counts[k,eps] for k in ("nontrivial first bound","nontrivial second bound") for eps in (0,1)):
        raise ArithmeticError(("missing positive controls",counts))
    print(dict(counts))


if __name__ == "__main__":
    verify()
