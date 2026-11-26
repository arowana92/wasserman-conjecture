"""Check the part bound, the leading-block tests and two large prime certificates."""
from collections import Counter
from itertools import combinations
from math import comb, gcd

from check_scaled_row_blocks import certify_prime
from search_endpoint import atoms, candidate, smallest_factors


def rows(n):
    for a in range(1,n//3+1):
        for b in range(a,(n-a)//2+1):
            yield (a,b,n-a-b), comb(n,a)*comb(n-a,b)


def layers(n, spf, predecessor_only=False):
    bases = atoms(n-1,spf)
    if not predecessor_only:
        bases += atoms(n,spf)
    for p,_ in bases:
        power = p
        while power <= n:
            if not predecessor_only or (n-1) % power == 0:
                yield p,power
            power *= p


def quantitative_controls():
    spf = smallest_factors(50)
    count = partial = 0
    for n in range(3,51):
        values = list(rows(n))
        for a, val_a in values:
            for b, val_b in values:
                G = gcd(n*(n-1),gcd(val_a,val_b))
                if a[2] == 1:
                    continue
                z = min(i for i in a if i >= 2)
                for _,q in layers(n,spf,True):
                    assert G*z >= q
                    partial += 1 < gcd(q,G) < q
                    count += 1
    assert partial
    print(f"{count} exact gcd bounds, including {partial} partially retained powers.")


def row_controls(limit):
    spf = smallest_factors(limit)
    counts = Counter()
    for n in range(3,limit+1):
        for row, value in rows(n):
            common = gcd(gcd(*row[:2]),row[2])
            scales = [u for u in range(3,common+1) if common%u == 0 and
                      gcd(u,n//u) == 1 and n//u >= 3]
            if not scales:
                continue
            avoiding = [(p,P) for p,P in layers(n,spf) if value%p]
            predecessor = [(q,Q) for q,Q in layers(n,spf,True) if value%q]
            for q,Q in predecessor:
                for u in scales:
                    assert min(row) > Q and n//u >= 2*Q+1
                    small = [a for a in row if a <= 2*Q]
                    assert not small or small == [Q+1] and (Q+1)%u == 0
                    counts["scale and six-part bounds"] += 1
                    counts["shifted boundary attained"] += bool(small)
            for p,P in avoiding:
                d,h = divmod(n,P)
                lead = [a//P for a in row]
                rem = [a%P for a in row]
                x = [d*a-j*n for a,j in zip(row,lead)]
                assert sum(lead) == d and sum(rem) == h and max(map(abs,x)) <= d*h
                if h and all(lead):
                    for q,Q in predecessor:
                        if d*h < Q:
                            delta = [a%Q for a in row]
                            assert sum(delta) == 1
                            assert all(d*(a-f) == j*(n-1) for a,f,j in zip(row,delta,lead))
                            counts["exact shifted blocks"] += 1
                for u in scales:
                    for _,Q in layers(n,spf,True):
                        # When the opposite row takes Q, u is its normalized sum.
                        if 0 < d*h <= 2*Q < u:
                            assert all(v == 0 for v in x)
                            counts["exact unshifted blocks"] += 1
            for (p,P),(q,Q) in combinations(avoiding,2):
                if p == q:
                    continue
                d,h = divmod(n,P)
                e,k = divmod(n,Q)
                if gcd(d,e) == 1 and d*e*(h+k) < n and d*e*h*k < n:
                    assert row[0]+row[1] <= min(h,k)
                    counts["two bases in one row"] += 1
    # Two constructions exercise the branches that first occur above the small scan.
    n, u, Q, P = 903, 129, 41, 128
    row = (u,2*u,4*u)
    d,h = divmod(n,P)
    assert (n-1)%Q == 0 and gcd(u,n//u) == 1 and 0 < d*h <= 2*Q < u
    assert comb(n,row[0])*comb(n-row[0],row[1]) % 2
    assert all(d*a == (a//P)*n for a in row)
    counts["exact unshifted blocks"] += 1
    P,Q = 2**58,11**16
    n = 11*P+55
    row = (22,33,n-55)
    value = comb(n,22)*comb(n-22,33)
    d,h = divmod(n,P)
    e,k = divmod(n,Q)
    assert n%11 == 0 and gcd(11,n//11) == 1 and value%2 and value%11
    assert gcd(d,e) == 1 and d*e*max(h*k,h+k) < n
    assert sum(row[:2]) <= min(h,k)
    counts["two bases in one row"] += 1
    assert all(counts[key] for key in ("scale and six-part bounds","shifted boundary attained",
                                       "exact unshifted blocks","exact shifted blocks",
                                       "two bases in one row"))
    print(dict(counts))


def second_interval():
    tests = 0
    for Q,R in combinations(range(2,55),2):
        for d in range(1,Q):
            if R-Q == d:
                continue
            allowed = [{w*m+d*delta for m in (-1,0,1) for delta in (0,1)}
                       for w in (Q,R)]
            assert set.intersection(*allowed) == {0,d}
            tests += 1
    print(f"{tests} two-power residue intersections.")


def examples():
    spf = smallest_factors(5636)
    for n in (5635,5636):
        d,h = divmod(n,512)
        e,k = divmod(n,625)
        assert n*(n-1)%10 == 0 and gcd(d,e) == 1
        assert min(d,e,h,k) >= 3 and d*e*max(h*k,h+k) < n
        assert max(gcd(d,n),gcd(e,n)) <= 2
        assert candidate(n,spf) is None
    e = 11*13*19
    for d,P,Q in (
        (2,3**80,108817983047835986093107686275454247),
        (127,2**256,5412450237894181046016723476830419859727020033384827332387109404041963682257),
    ):
        for prime in (3,5,7,11,13,17,19,Q):
            certify_prime(prime)
        n = e*Q+1
        h = n-d*P
        assert P <= 2*e*h <= 2*P and h*h > n and n%1785 == 0
        assert n*(n-1)%(3 if d == 2 else 2) == 0
        if d == 2:
            assert h < Q
        else:
            assert d*h < Q and n%d == 2 and (d*h)**5 > 27*(n-1)**2
        print(f"Certified {n} with leading quotient {d} and offset {h}.")


if __name__ == "__main__":
    quantitative_controls()
    row_controls(250)
    second_interval()
    examples()
