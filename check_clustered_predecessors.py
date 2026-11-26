"""Independent binomial, CRT and adjacently shared certificate checks."""

from collections import Counter
from fractions import Fraction
from itertools import combinations, product
from math import comb, gcd, prod

from flint import fmpz

from check_cluster_certificates import PRIMES, PRODUCT, LIMIT, universal_prime
from check_cofactor_digit_gap import atoms
from check_minimality import digit_sum
from check_scaled_row_blocks import certify_prime
from check_three_factors import valuation


def binomial_avoids(n, a, p):
    power = p
    while power <= n:
        if n//power != a//power + (n-a)//power:
            return False
        power *= p
    return True


def exact_controls():
    counts = Counter()
    prime_list = [p for p, w in atoms(180) if w == p]
    for n in range(3, 180):
        rows = [(a, b, n-a-b) for a in range(1, n//3+1) for b in range(a, (n-a)//2+1)]
        coeffs = [comb(n, a)*comb(n-a, b) for a, b, _ in rows]
        u = universal_prime(n)
        common = gcd(*coeffs) if len(coeffs)>1 else coeffs[0]
        assert (common > 1) == bool(u)
        if u:
            assert common % u == 0
        for p, q in combinations([v for v in prime_list if (n-1)%v == 0], 2):
            m = n-1
            joint = any(comb(m, a)%p and comb(m, a)%q for a in range(p*q, m, p*q))
            for (a, b, c), k in zip(rows, coeffs):
                if k%p and k%q:
                    assert joint
                    z = next(x for x in (a,b,c) if x%(p*q) == 0)
                    assert 0 < z < m and gcd(comb(m,z),p*q) == 1
                    counts['successful transfers'] += 1
                counts['successor comparisons'] += 1
    return counts


def first_shared_minimum():
    previous = set()
    for n in range(3,53):
        coeffs = [comb(n,a)*comb(n-a,b) for a in range(1,n//3+1)
                  for b in range(a,(n-a)//2+1)]
        edges = set()
        if gcd(*coeffs) == 1:
            primes = [v for v in range(2,n+1) if fmpz(v).is_prime()]
            masks = {p: sum(1 << i for i,k in enumerate(coeffs) if k%p) for p in primes}
            edges = {(p,q) for p,q in combinations(primes,2) if not (masks[p]&masks[q])}
        common = previous & edges
        vertices = sorted({p for pair in common for p in pair})
        lists = [group for group in combinations(vertices,4)
                 if all(pair in common for pair in combinations(group,2))]
        if lists:
            assert n == 52 and (2,3,5,17) in lists
            return n-1, n, lists
        previous = edges
    raise AssertionError('shared minimum missing')


def general_binomial_checks():
    universe = atoms(100)
    counts = Counter()
    for (p, W), (q, Z) in combinations(sorted(universe, key=lambda v: v[1]), 2):
        if p == q:
            continue
        for m in range(2, 45):
            n = W*Z*m
            # Enumerate all normalized positive binomial rows independently of the bounds.
            joint = any(binomial_avoids(n, W*Z*f, p) and binomial_avoids(n, W*Z*f, q)
                        for f in range(1, m//2+1))
            if joint:
                counts['jointly avoidable weights and pairs'] += 1
            for K in (Z, 2*Z, 3*Z):
                center = W*Z*K
                lw = 1
                while lw*p <= center:
                    lw *= p
                lz = q
                while lz <= center:
                    lz *= q
                for lower in (lw, lw//p):
                    if lower < p:
                        continue
                    D, H = divmod(n, lower)
                    E = n//lz
                    first = D == E and (-W%Z)*(m-K*D) <= Z
                    second = E >= D-1 and (-W%Z)*(H//(W*Z)) + (n//Z)%Z < Z
                    if first or second:
                        assert not joint, (n, W, Z, K, lower, lz)
                    counts['first criterion'] += first
                    counts['second criterion'] += second
    return counts


def cluster_checks(limit):
    ascending = sorted(atoms(limit), key=lambda v: v[1])
    counts = Counter()
    for i, (p, P) in enumerate(ascending):
        for j in range(i+1, len(ascending)):
            q, Q = ascending[j]
            if 2*(Q-P)**2 > P:
                break
            if q == p:
                continue
            for r, R in ascending[j+1:]:
                h = R-P
                if 2*h*h > P:
                    break
                if r in (p, q):
                    continue
                ps, ws = (p, q, r), (P, Q, R)
                M = P*Q*R
                for d in (1, 2):
                    if 2*d*h*h > P or (d == 2 and M%2 == 0):
                        continue
                    n = d*M+1
                    inverses = [pow(M//w, -1, w) for w in ws]
                    for extra in product(range(d), repeat=3):
                        k = [a+j*w for a, j, w in zip(inverses, extra, ws)]
                        row = [(M//w)*b for w, b in zip(ws, k)]
                        if sum(row) != n:
                            continue
                        counts['canonical CRT rows'] += 1
                        if valuation(n,row,p) == 0:
                            floors = [a//P**3 for a in row]
                            assert sum(k) == sum(b*w for b,w in zip(floors,ws))
                            assert d == 1 and gcd(*k) == 1 and M%2 == 0
                            counts['rounded survivors'] += 1
                if R-Q <= Q-P:
                    if all(b == w for b,w in zip(ps,ws)) and R-Q >= 3:
                        f = R//(R-Q)
                        row = [Q*R*f]*2 + [Q*R*(P-2*f)]
                        assert min(row)>0 and not valuation(M,row,q) and not valuation(M,row,r)
                        counts['opposite-gap compatible rows'] += 1
                    continue
                counts['increasing-gap triples'] += 1
                for d in range(1,P//(2*h*h)+1):
                    n = d*M
                    for (a,W), (b,Z) in combinations(zip(ps,ws),2):
                        for f in range(1,n//(W*Z)//2+1):
                            part = f*W*Z
                            assert not (binomial_avoids(n,part,a) and binomial_avoids(n,part,b))
                            counts['normalized binomial rows'] += 1
    return counts


def fourth_prime_checks():
    counts = Counter()
    powers = atoms(120)
    for n in range(6,160):
        for p, W in powers:
            if n%W:
                continue
            for r, L in powers:
                if r == p or L > n:
                    continue
                m, a = divmod(n, L)
                if a+1 >= L or m*(a+1) >= W:
                    continue
                g = gcd(m,a)
                for shift in (0, 1):
                    weight = n+shift
                    counts['tested classes'] += 1
                    # The W-residues give parts W*x+delta in both weights.
                    for x in range(n//W+1):
                        for y in range(n//W-x+1):
                            z = n//W-x-y
                            row = (W*x, W*y, W*z+shift)
                            if not min(row) or valuation(weight,row,p) or valuation(weight,row,r):
                                continue
                            assert g > 1
                            assert all(part % (n//g) == 0 for part in (W*x,W*y,W*z))
                            counts['compatible rows at noncoprime quotients'] += 1
    return counts


def close_product_checks():
    Y, H = 40000, 10
    data = []
    for P in range(Y,Y+100):
        for s in range(1,H):
            for t in range(1,H-s+1):
                Q, R = P+s, P+s+t
                S = P+Q+R
                B = s*s+s*t+t*t
                C = -(2*s+t)*(s-t)*(s+2*t)
                M = P*Q*R
                assert 27*M == S**3-3*S*B+C
                data.append((M,S,B,C))
    data.sort()
    counts = Counter()
    for i, (M,S,B,C) in enumerate(data):
        j = i
        while j < len(data) and data[j][0]-M <= Y:
            other, total, norm, cubic = data[j]
            assert total == S
            assert 27*(other-M) == 3*S*(B-norm)+cubic-C
            counts['equal norms' if B == norm else 'distinct norms'] += 1
            j += 1
    assert all(counts.values()) and len(counts) == 2
    return counts


def separated_examples():
    for B,X,Y,a,b,first in ((2**32,2**24,2**32,43,15,True),
                             (2**20,2**20,2**90,7,133,False)):
        certify_prime(X+a)
        certify_prime(Y+b)
        C = max(B,X,Y)
        v, w = max(((a,X),(b,Y)), key=lambda f: Fraction(*f))
        if first:
            assert 8*C*v*v <= w*w
        else:
            assert B <= min(X,Y) and 8*a*b*v <= w


def check_minima():
    # Factor the full endpoint, which contains every individually universal prime.
    ds = list(range(1,501)) + [951125, LIMIT]
    for d in ds:
        n = d*PRODUCT
        factorizations = [fmpz(k).factor() for k in (n-1,n,n+1)]
        for i in (0,1):
            primes = {int(p) for factors in factorizations[i:i+2] for p, _ in factors}
            found = {p for p in primes if digit_sum(n+i,p) <= 2}
            witness = universal_prime(n+i)
            assert witness in found if found else witness == 0
    for p in PRIMES:
        certify_prime(p)
    for d in (1, LIMIT):
        n = d*PRODUCT
        near = n if n%2 else n-1
        while not fmpz(near).is_prime():
            near -= 2
        certify_prime(near)
        assert n+1-near <= PRIMES[0]
        print(f'd={d}, weight {n}, fourth prime {near}, minima '
              f'{1 if universal_prime(n) else 4} and {1 if universal_prime(n+1) else 4}')
    print(f'Factored endpoint checks of minimum size: {2*len(ds)}')


if __name__ == '__main__':
    print('Exact coefficients:', exact_controls())
    print('First shared minimum:', first_shared_minimum())
    print('Binomial lemma:', general_binomial_checks())
    print('Prime-power triples:', cluster_checks(1200))
    print('Fourth-prime test:', fourth_prime_checks())
    print('Close integer products:', close_product_checks())
    separated_examples()
    check_minima()
