"""Verify the minor inequalities, the exact-row test and lowest-layer examples."""
from collections import Counter
from itertools import product
from math import comb, gcd

from flint import fmpz

from check_central_cuts import legendre
from check_scaled_row_blocks import certify_prime


def sorted_rows(n):
    for a in range(1,n//3+1):
        for b in range(a,(n-a)//2+1):
            yield (a,b,n-a-b),comb(n,a)*comb(n-a,b)


def minors(row,r):
    return [r[j]*row[k]-r[k]*row[j] for j,k in ((1,2),(2,0),(0,1))]


def coefficient_controls():
    counts = Counter()
    for n in range(3,39):
        rows = list(sorted_rows(n))
        for r,u in (((1,2,4),(2,3,5)),((1,2,3),(1,1,1)),((2,3,1),(4,2,1)),
                    ((1,1,0),(1,2,0)),((1,1,1),(2,0,0))):
            left = [(a,v,minors(a,r)) for a,v in rows]
            right = [(b,v,minors(b,u)) for b,v in rows]
            for (row_a,v,ma),(row_b,w,mb) in product(left,right):
                common = gcd(v,w)
                D,J = gcd(n,common),gcd(n-1,common)
                X,Y = max(map(abs,ma)),max(map(abs,mb))
                if all(ma) and all(mb):
                    assert n**3*(n-1) <= D**3*J*(X*Y)**3
                    counts["six nonzero minors"] += 1
                if X and all(mb):
                    assert 2*n*n*(n-1) < D**3*J*X*X*Y**3
                    counts["five nonzero minors"] += 1
                    counts["exactly five"] += 0 in ma
                if Y:
                    h = sum(row_a[:2])
                    assert 2*n*n*(n-1) < D**3*J*h**3*Y**2
                    counts["two small parts and a nonexact companion"] += 1
                if not any(ma):
                    s = sum(r)
                    assert n%s == 0
                    Gamma = gcd(r[1],r[2])*gcd(r[2],r[0])*gcd(r[0],r[1])
                    # The part from M unshared by the two coefficients must
                    # divide Gamma times the row-B zero-or-one products.
                    assert all(Gamma*J*b*(b-1) % (n-1) == 0 for b in row_b)
                    counts["exact target rows"] += 1
    assert counts["exactly five"]
    print(dict(counts))


def square_controls():
    counts = Counter()
    for t in (3,5,7):
        for T in range(2,1300):
            residues = [j for j in range(T) if j*(j-t)%T == 0]
            values = {r+h*T for r in residues for h in (-1,0,1)}
            for a,b in product(values,repeat=2):
                c = t-a-b
                if c not in values:
                    continue
                z = (2*a-t,2*b-t,2*c-t)
                Z = max(map(abs,z))
                assert sum(z) == -t and all((v*v-t*t)%(4*T) == 0 for v in z)
                if T >= 16*t*t and 27*t*Z**3 < 8*T*T:
                    assert -t in z
                    counts["strict window"] += 1
                    counts["one unshifted only"] += z.count(-t) == 1
                if all(v*v > t*t for v in z):
                    counts["all three positive square quotients"] += 1
    # A small nonconstant quadratic inside the strict test.
    # Choose T = b*(b-t) so that the triple of shifted variables is (0,b,t-b).
    t,b = 7,10000
    T = b*(b-t)
    z = (-t,2*b-t,t-2*b)
    assert T >= 16*t*t and 27*t*max(map(abs,z))**3 < 8*T*T
    assert all((v*v-t*t)%(4*T) == 0 for v in z) and z.count(-t) == 1
    counts["one unshifted only"] += 1
    assert counts["strict window"] and counts["all three positive square quotients"]
    print(dict(counts))


def progression_controls():
    levels = 0
    for h in (*range(20),10**6,10**18):
        v = 11972031+75150075*h
        u = 2*v-1
        n,S,T = u*v,v-1,2*v+1
        m,l = S//35,T//143
        A = (u*(23*m+1),u*5*m,u*7*m)
        B = (v*(119*l-2),v*11*l,v*13*l)
        assert u-2*v == -1 and S*T == n-1 and gcd(u,v) == gcd(S,T) == 1
        assert S == 35*m and T == 143*l and m%35 == 3 and l%143 == 131
        for scale,row,shift in ((u,A,0),(v,B,0),(S,A,1),(T,B,1)):
            assert sum(row) == n
            for p,e in fmpz(scale).factor():
                w = int(p)**int(e)
                assert sum(a%w for a in row) == shift == n%w
                levels += 1
        for row in (A,B):
            assert legendre(n,2) > sum(legendre(a,2) for a in row)
    print(f"{levels} full endpoint layers in 22 progression members.")


def whole_weight_controls():
    P = 5**168
    q = 42211883872453072075928572096340406768795975134335019729909
    certify_prime(q)
    certify_prime(5)
    Q = q*q
    lower,upper = 200_000_000,300_000_000
    lo,hi = (q*(3*q+m) for m in (lower,upper))
    h0,h1 = lo-2*P,hi-2*P
    k0,k1 = lower*q,upper*q
    assert 0 < h0 <= h1 < P and 0 < k0 <= k1 < Q
    assert min(h0,k0)**7*1000 > hi**4 and max(h1,k1)**7 < lo**4
    assert max(h1,k1)**7 <= 6**6*(lo-1)**4 and 6*(h1+k1) < lo
    assert max((2*h1)**2*(3*k1)**3,(2*h1)**3*(3*k1)**2) <= 2*lo*lo*(lo-1)
    assert lo-1 >= 64 and 54*(4*h1+2)**3 < 8*(lo-1)**2
    classes = [r for r in range(5) if q*(3*q+r)%5 in (0,1)]
    count = sum((upper-r)//5-(lower-1-r)//5 for r in classes)
    assert classes == [2,3] and count == 40_000_000
    n = q*(3*q+256640162)
    assert len(str(n)) == 118 and (n-1)%(5*7*11*13) == 0 and n%17 == 0
    W,m = 5**16,132147841523
    n = q*(3*q+m)+1
    h,k = n-2*P,n-3*Q
    assert 0 < m < W and (n-1)%(q*W) == 0 and k <= q*W
    assert min(h,k)**7 > 6**6*(n-1)**4 and 28282*min(h,k)**5 > n**3
    assert 6*(h+k) < n and 0 < h < P and 0 < k < Q
    assert max((2*h)**2*(3*k)**3,(2*h)**3*(3*k)**2) <= 2*n*n*(n-1)
    assert n-1 >= 64 and 54*(4*h+2)**3 < 8*(n-1)**2
    for n in range(3,251,2):
        z = n//2
        for a in range(1,(n-z)//2+1):
            assert comb(n,z)*comb(n-z,a)%2 == 0
    print(f"{count} weights in the two- and three-block example.")


if __name__ == "__main__":
    coefficient_controls()
    square_controls()
    progression_controls()
    whole_weight_controls()
