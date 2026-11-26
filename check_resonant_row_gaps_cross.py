"""Second enumeration of resonant gap shapes, using coprime products and digit sums."""

from argparse import ArgumentParser
from collections import Counter
from functools import lru_cache
from itertools import combinations
from math import gcd, isqrt, prod


@lru_cache(maxsize=8192)
def factors_of(n):
    return sorted({k for i in range(1,isqrt(n)+1) if n % i == 0 for k in (i,n//i)})


def splits(n):
    return ((k,n//k) for k in factors_of(n) if gcd(k,n//k) == 1)


def bin_systems(t):
    for thirds in factors_of(t):
        for mids in factors_of(2*t-2):
            if gcd(thirds,mids) > 1:
                continue
            for r,z in splits(thirds):
                for s,w in splits(mids):
                    if (s == 1 or r == 1) and (w == 1 or z == 1):
                        continue
                    left, right = s*r, w*z
                    # b/left and c/right are reciprocal; only a factor 2 can cancel.
                    for alpha,beta in ((1,1),(1,2),(2,1)):
                        if left % alpha or right % beta:
                            continue
                        b,c = alpha*(right//beta), beta*(left//alpha)
                        if b < 2*r+s or c < 2*z+w:
                            continue
                        for x,y in splits(t):
                            if x % r == y % z == 0 and (t+1) % gcd(x*c,y*b) == 0:
                                yield x,y,r,z,s,w,b,c


def rows(total, divisor, stride):
    for third in range(divisor,total,divisor):
        remainder = total-third
        if remainder % (2*stride):
            continue
        unit = remainder//2
        if gcd(unit,divisor) == gcd(third,stride) == 1:
            yield unit,third


def digits_pass(parts, p, bound):
    parts = [k % bound for k in parts]
    while bound > 1:
        if sum(k % p for k in parts) >= p:
            return False
        parts = [k//p for k in parts]
        bound //= p
    return True


def numeric(x,y,b,c,m,h,n,j,S,T):
    U,V = c*T-y,b*S+x
    return U,V,(U*(m*S+x),U*m*S,U*h*S),(V*n*T,V*(n*T-y),V*j*T)


@lru_cache(maxsize=8192)
def possible_leading_digits(params, p, S, T):
    # Lift the two residues to a positive integral solution. Exact division of
    # its actual parts gives the digit slopes, which depend only on S,T mod p.
    x,y,b,c,*_ = params
    g = gcd(x*c,y*b)
    a,v,target = x*c//g,y*b//g,(x*y+1)//g
    correction, rem = divmod(a*T-v*S-target,p)
    assert rem == 0
    k = (correction*pow(v,-1,a)) % a + a
    S += p*k
    T = (v*S+target)//a
    U,V,A,C = numeric(*params,S,T)
    for factor,row,exception in ((U,A,-1),(V,C,-1),(S,A,0),(T,C,1)):
        if factor % p:
            continue
        slopes = []
        for i,value in enumerate(row):
            quotient,remainder = divmod(value-int(i == exception),factor)
            assert remainder == 0
            slopes.append(quotient % p)
        allowed = set(range(1,p))
        for unit in tuple(allowed):
            if sum(unit*k % p for k in slopes) >= p:
                allowed.remove(unit)
        if not allowed:
            return False
    return True


def test_modulus(params,bins,p,q):
    x,y,b,c,*_ = params
    r,z,s,w = bins
    common = gcd(x*c,y*b)
    a, v, target = x*c//common, y*b//common, (x*y+1)//common
    # Solve the divided relation directly for one bin modulo q.
    solve_for_second = gcd(a,p) == 1
    inverse = pow(a if solve_for_second else v,-1,q)
    for residue in range(q):
        if solve_for_second:
            S,T = residue, (inverse*(target+v*residue)) % q
        else:
            S,T = (inverse*(a*residue-target)) % q, residue
        if not possible_leading_digits(params,p,S % p,T % p):
            continue
        U,V,A,C = numeric(*params,S,T)
        if sum(j % p == 0 for j in (*bins,S,T)) > 1 or U % p == V % p == 0:
            continue
        if U*S*r*s % p == 0 and not digits_pass(A,p,q):
            continue
        if V*T*w*z % p == 0 and not digits_pass(C,p,q):
            continue
        return False
    return True


def finite_bins(params,first):
    x,y,b,c,m,_,n,_ = params
    if first:
        for S in range(1,(b*x-1)//m+1):
            T,rem = divmod(x*y+1+y*b*S,x*c)
            if rem == 0:
                yield S,T
    else:
        for T in range(1,((c+1)*y-1)//n+1):
            S,rem = divmod(x*c*T-x*y-1,y*b)
            if rem == 0 and S > 0:
                yield S,T


def audit(limit, first_t=2, check_remaining=None):
    counts = Counter()
    maximum = 0
    for t in range(first_t,limit+1):
        for x,y,r,z,s,w,b,c in bin_systems(t):
            for m,h in rows(b,s,r):
                for n,j in rows(c,w,z):
                    counts['shapes'] += 1
                    params = x,y,b,c,m,h,n,j
                    bins = r,z,s,w
                    for first in (False,True):
                        if (w*z if first else s*r) % 2 == 0:
                            continue
                        if not first and j % (2*n) == 0:
                            # A fixed list suffices for the reported interval.
                            if not any(test_modulus(params,bins,p,q)
                                       for p,q in ((3,81),(5,25),(19,361),(2,8192))):
                                if check_remaining is None:
                                    raise ArithmeticError(('unresolved progression',params,bins))
                                check_remaining(params,bins)
                                counts['checked to exponent bound'] += 1
                            else:
                                counts['unbounded'] += 1
                            continue
                        for S,T in finite_bins(params,first):
                            U,V,A,C = numeric(*params,S,T)
                            if min(*A,*C) <= 0:
                                continue
                            assert U*V == sum(A) == sum(C) == prod(bins)*S*T+1
                            counts['positive'] += 1
                            maximum = max(maximum,U*V)
                            if gcd(U,V) > 1 or any(gcd(i,j) > 1 for i,j in combinations((*bins,S,T),2)):
                                continue
                            counts['coprime'] += 1
                            if digits_pass(A if first else C,2,1 << (U*V).bit_length()):
                                raise ArithmeticError(('requires further tests',U*V,A,C))
    print('Independent counts',dict(counts),'largest finite weight',maximum)


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('limit', type=int, nargs='?', default=2660)
    audit(parser.parse_args().limit)
