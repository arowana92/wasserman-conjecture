"""Search the forced near-equality bin triples by divisors of a cubic discriminant."""

from argparse import ArgumentParser
from collections import Counter
from itertools import combinations
from math import gcd, isqrt, prod

from flint import fmpz
from check_three_factors import valuation


def factors(n):
    values = [(int(p), e) for p, e in fmpz(n).factor()]
    if prod(p**e for p, e in values) != n or any(not fmpz(p).is_prime() for p, e in values):
        raise ArithmeticError(n)
    return values


def divisors(factorization):
    values = [1]
    for p, e in factorization:
        values = [d*p**j for d in values for j in range(e+1)]
    return values


def check(h, k, counts):
    s, b = h+k, h*k
    disc = (b*(k-h))**2 + 4*s**3 - 18*s*b - 27
    if not disc:
        raise ArithmeticError((h, k))
    counts['gap pairs'] += 1
    for u in divisors(factors(abs(disc))):
        if u <= b or u % 4 != 3:
            continue
        root = isqrt(s*s+3*(u-b))
        if root**2 != s*s+3*(u-b) or (root-s) % 3:
            continue
        d = (root-s)//3
        e, f = d+h, d+k
        if d <= 1 or d % 2 == 0 or gcd(d,e)*gcd(d,f)*gcd(e,f) != 1:
            continue
        t = d*e*f
        if (t-1) % u or (t-1)//u < 2 or gcd(t,u-1) != 1:
            continue
        counts['bin triples'] += 1
        n = (u-1)*t+1
        v = n//u
        assert u == d*e+d*f+e*f and n == u*v
        print('BINS', d, e, f, u, n)
        if gcd(u,v) != 1:
            continue
        counts['unitary weight'] += 1
        if any(valuation(n,(u,n-u),p) for p,_ in factors(u-1)+factors(u)):
            continue
        counts['anchor binomial'] += 1
        row = tuple(v*y for y in (d*e,d*f,e*f))
        if any(valuation(n,row,p) for p,_ in factors(t)+factors(v)):
            continue
        counts['full companion'] += 1
        ell = (t-1)//u
        step = u*(u-1)
        for j in range(1,ell//2+1):
            row_a = (u,step*j,step*(ell-j))
            if all(valuation(n,row_a,p) == 0 for p,_ in factors(u)+factors(u-1)):
                raise ArithmeticError((n, row_a, row))


def audit_bins(path):
    counts = Counter()
    largest = (0, 0)
    with open(path) as stream:
        for line in stream:
            if not line.startswith('BINS '):
                continue
            d,e,f,u,n = map(int,line.split()[1:])
            t = d*e*f
            ell, rem = divmod(t-1,u)
            assert not rem and u == d*e+d*f+e*f and n == (u-1)*t+1
            m = d+e+f-9*ell
            assert (m-3) % 8 == (u-3) % 8 == 0
            assert all((v-ell) % 4 == 1 for v in (d,e,f))
            y = (d*e,d*f,e*f)
            H = sum(x*(z-w)**2 for x,z,w in ((d,e,f),(e,f,d),(f,d,e)))
            assert m*u == H+9 and H >= 3*u-9
            deficit = u**3-27*t**2
            weighted = sum(x*(z-w)**2 for x,z,w in ((y[0],y[1],y[2]),
                                                       (y[1],y[2],y[0]),(y[2],y[0],y[1])))
            assert 2*deficit == u*sum((r-s)**2 for r,s in combinations(y,2))+6*weighted
            row = tuple(n//u*v for v in y)
            carries = valuation(n,(u,n-u),2), valuation(n,row,2)
            # Independent binary test for the first binomial and the explicit companion.
            assert bool(u & (n-u)) == bool(carries[0])
            assert bool((row[0]&row[1]) | (row[1]&row[2]) | (row[2]&row[0])) == bool(carries[1])
            counts['triples'] += 1
            counts['both even'] += min(carries)>0
            largest = max(largest,(n,u))
    print(dict(counts), 'largest weight and part', largest)


def check_prime_power_blocks():
    cases = positives = 0
    for p in (2,3,5,7,11,13):
        U = p
        while U <= 130:
            for x in (1,p,p*p):
                a = U*x
                for k in (1,2):
                    if (a-1) % k:
                        continue
                    b = (a-1)//k
                    for l in range(2,U):
                        for j in range(1,l//2+1):
                            row = (a, U*b*j, U*b*(l-j))
                            unit = valuation(sum(row),row,p)==0
                            cases += 1
                            if unit:
                                positives += 1
                                assert k == 2 and j % 2 and (l-j) % 2
            U *= p
    print('prime-power block tests', cases, 'carry-free cases', positives)


def check_lattice_residues():
    counts = Counter()
    # Fix two odd bins. The equation T-ell*U=1 determines the third.
    for coefficient in (1,2):
        for d in range(1,256,2):
            for e in range(1,256,2):
                den_slope = coefficient*e+d
                for ell in range(1,(d*e-1)//den_slope+1):
                    f, leftover = divmod(1+ell*d*e,d*e-ell*den_slope)
                    if leftover or f%2==0:
                        continue
                    u = coefficient*e*f+d*f+d*e
                    assert d*e*f-ell*u==1
                    if coefficient==1:
                        # With odd bins U is odd, so the equation forces even ell.
                        assert ell%2==0 and u%8==3
                        assert all((b-ell)%4==1 for b in (d,e,f))
                    else:
                        assert ell%2==0 and u%8 in (2,4)
                    counts[coefficient]+=1
    print('lattice residue examples with determinant one and two',counts[1],counts[2])


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('gap', type=int, nargs='?')
    parser.add_argument('--audit')
    parser.add_argument('--blocks', action='store_true')
    args = parser.parse_args()
    if args.audit:
        audit_bins(args.audit)
    if args.blocks:
        check_prime_power_blocks()
        check_lattice_residues()
    if args.gap is None:
        if not args.audit and not args.blocks:
            parser.error('supply a gap, --audit, or --blocks')
        return
    counts = Counter()
    for h, k in combinations(range(2,args.gap+1,2),2):
        check(h,k,counts)
    print(dict(counts))


if __name__ == '__main__':
    main()
