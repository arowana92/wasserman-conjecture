"""Independent enumeration of successor assignments by weight and checks using exact coefficients."""

from itertools import combinations
from math import comb, gcd, isqrt, prod

from check_successor_powers import first_rows, square
from search_endpoint import atoms, smallest_factors


def products(values):
    return [prod(s) for size in range(len(values)+1) for s in combinations(values, size)]


def compare_assignments(limit, minimum):
    cutoff = isqrt(limit**5//27)
    # Factor every Q by a sieve. For N we need only its full powers with base < d.
    spf = smallest_factors(cutoff+limit)
    primes = [r for r in range(2, limit) if spf[r] == r]
    direct = set()
    for P in range(2, cutoff+1):
        pf = atoms(P, spf)
        if len(pf) != 1:
            continue
        p, = (r for r, _ in pf)
        for d in range(4, min(P, limit+1)):
            Q = P+d
            n = P*Q+1
            if n <= minimum or 27*Q*Q >= d**5:
                continue
            factors = atoms(Q, spf)
            if max(w for _, w in factors) > d:
                continue
            available = []
            for r in primes:
                if n % r:
                    continue
                w = r
                while n % (r*w) == 0:
                    w *= r
                if w < d:
                    available.append(w)
            for u in products(available):
                if not 3 <= u < d:
                    continue
                eligible = [w for r,w in factors if r != p and 2*w+1 <= u]
                for T in products(eligible):
                    if 27*T*T > u**3:
                        continue
                    B = Q//T
                    for y in range(u,d,u):
                        e = 1
                        for _, w in factors:
                            if T % w == 0 or (y//u) % w == 0:
                                continue
                            if (P*y-1) % w:
                                break
                            e *= w
                        else:
                            if e <= d:
                                direct.add((P,d,y,u,B,e))
    forward = set()
    square(limit, minimum, forward)
    if direct != forward or not direct:
        raise ArithmeticError(("assignments", len(direct), len(forward), list(direct ^ forward)[:5]))
    print(f"Weight-first and progression enumerations agree on {len(direct)} assignments.")


def compare_exact(maximum):
    spf = smallest_factors(2*maximum**2+1)
    total = scaled = exceptional = 0
    for P in range(2, maximum+1):
        pf = atoms(P, spf)
        if len(pf) != 1:
            continue
        p = pf[0][0]
        for d in range(2,P):
            Q = P+d
            n = P*Q+1
            for t in range(2,d+1):
                for y in range(1,t):
                    x = t-y
                    for row in ((P*y,P*x,P*(Q-t)+1), (P*y,P*x+1,P*(Q-t))):
                        a,b,_ = row
                        coefficient = comb(n,a)*comb(n-a,b)
                        if coefficient % p == 0:
                            continue
                        total += 1
                        u = prod(w for r,w in atoms(n,spf) if coefficient % r)
                        B = prod(w for r,w in atoms(Q,spf) if coefficient % r)
                        e = gcd(B,P*y-1)
                        if y % u or (y//u) % (B//e) or (d*y+1) % e or e > d:
                            raise ArithmeticError(("exact divisibilities",n,row,u,B,e))
                        if row not in first_rows(P,Q,d,y,u,e):
                            raise ArithmeticError(("missing exact row",n,row,u,e))
                        scaled += u >= 3
                        exceptional += e > 1
    if not (scaled and exceptional):
        raise ArithmeticError((total, scaled, exceptional))
    print(f"{total} exact p-free rows, including {scaled} with U >= 3 and {exceptional} with e > 1.")


if __name__ == "__main__":
    compare_assignments(60, 10000)
    compare_exact(42)
