"""Bound all non-root shapes and test the companion reconstruction without factoring H."""

from itertools import combinations_with_replacement
from math import gcd, isqrt, lcm


def classify():
    seen = 0
    roots = set()
    upper = 0
    for a in range(11):
        for b in range(8):
            if not (a or b):
                continue
            for c in range(11-a):
                for e in range(8-b):
                    if (c,e) < (a,b) or (a+c,b+e) == (10,7):
                        continue
                    seen += 1
                    bound = gcd(a*a-7*a*b+10*b*b, c*c-7*c*e+10*e*e, a*e-b*c)
                    if not bound:
                        assert b and e and a % b == 0 and c*b == a*e
                        u = a//b
                        roots.add((7-u,u,b,e))
                        continue
                    available = lcm(abs(10*b-7*a), abs(10*e-7*c))
                    assert available
                    if bound >= 3:
                        upper = max(upper, (isqrt(bound**3*available**2//27)-7)//10)
    expected = {(t,7-t,i,j) for t in (2,5) for i in range(1,t//2+1) for j in range(i,t-i+1)}
    assert roots == expected and (seen,upper) == (1069,3981)
    assert upper*(10*upper+7)+1 < 10**9
    print(f"{seen} shapes, {len(roots)} root shapes, bounded cases end by P={upper}.")


def exact_parts(bound):
    counts = {}
    for t in (2,5):
        total = even = 0
        for P in range(1,bound+1):
            if t == 5 and P%2 == 0 and P & (P-1):
                continue
            T,W,H = 1+t*P, 1+(7-t)*P, 10*P+7
            # A normalized part is an idempotent after multiplying by T
            # at each prime-power factor. Test the single product directly.
            possible = [f for f in range(1,W-1) if (T*f)*(T*f-1) % H == 0]
            lookup = set(possible)
            for a,b in combinations_with_replacement(possible, 2):
                c = W-a-b
                if c < b or c not in lookup:
                    continue
                total += 1
                x,y,z = T*a,T*b,T*c
                assert t == 2 or P%2
                assert (x&y) or (y&z) or (x&z)
                even += 1
        counts[t] = total
        assert even == total
    assert counts == {2:58, 5:1}, counts
    print(f"Factor-free part enumeration through {bound}: {counts} rows.")


if __name__ == "__main__":
    classify()
    exact_parts(2500)
