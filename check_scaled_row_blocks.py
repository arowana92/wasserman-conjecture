"""Verify a digit obstruction for scaled rows and an explicit optimal four-prime set."""

from functools import lru_cache
from math import comb, gcd, prod

from flint import fmpz

from search_endpoint import atoms, smallest_factors


def block_fails(n, scale, p, w):
    target = n//scale
    m = p
    while m <= p*n:
        negative = (-scale)%m
        if negative*((target-1)//w) < 2*m-(n//w)%m:
            return True
        m *= p
    return False


def exact_controls():
    spf = smallest_factors(150)
    row_tests = nontrivial = rejected = allowed = 0
    for n in range(3,151):
        rows = [(a,b,n-a-b,comb(n,a)*comb(n-a,b)) for a in range(1,n//3+1)
                for b in range(a,(n-a)//2+1)]
        for scale in range(1,n//3+1):
            if n%scale:
                continue
            target = n//scale
            for p, full in atoms(n-1,spf):
                w = p
                while w <= full:
                    fails = block_fails(n,scale,p,w)
                    if target >= 2*w+1:
                        nontrivial += fails
                    for a,b,c,value in rows:
                        if a%scale or b%scale or c%scale:
                            continue
                        row_tests += 1
                        if fails:
                            assert value%p == 0, (n,scale,p,w,a,b)
                            rejected += 1
                        else:
                            allowed += value%p != 0
                    w *= p
    assert min(nontrivial,rejected,allowed) > 0
    print(f'{row_tests} exact coefficient tests, {rejected} certified divisibilities, '
          f'{allowed} prime avoidances and {nontrivial} obstructions beyond the size test.')


def square_root_controls():
    spf = smallest_factors(20000)
    first = companion = 0
    for k,s,t in ((3,1,1),(2,1,4),(2,2,1),(4,1,2),(3,1,3),(5,1,1)):
        coefficient = (k+1)*s
        for g in range(3,1000):
            if gcd(t,g) != 1:
                continue
            d = coefficient*g+t
            if k*s%t == 0:
                for r,w in atoms(g,spf):
                    if w*w > coefficient*g:
                        digit = -s*(g//w)*pow(t,-1,w)%w
                        assert digit + (k*digit)%w >= w
                        first += 1
            if s%t == 0:
                for r,w in atoms(d,spf):
                    if t*w*w >= d and 2*w+1 <= g:
                        assert (coefficient//t%w)*((g-1)//w) < w+1
                        companion += 1
    assert first and companion
    print(f'{first} first-row and {companion} companion square-root checks.')


@lru_cache(maxsize=None)
def certify_prime(n):
    # FLINT suggests a factorization. The following order and gcd checks
    # recursively certify primality over the integers.
    if n == 2:
        return
    powers = [(int(q),int(e)) for q,e in fmpz(n-1).factor()]
    assert all(q >= 2 and e >= 1 for q,e in powers)
    assert prod(q**e for q,e in powers) == n-1
    for q,_ in powers:
        certify_prime(q)
        a = 2
        while gcd(pow(a,(n-1)//q,n)-1,n) != 1 or pow(a,n-1,n) != 1:
            a += 1
            if a == n:
                raise ArithmeticError(('no primality witness',n,q))
    # Every prime divisor of n is 1 modulo every full power in n-1.


def four_primes():
    w = 2**127-1
    p,q = 80*w-1,282*w-1
    n = p*q
    r = n-94
    for number in (w,p,q,r):
        assert fmpz(number).is_prime()
        certify_prime(number)
    assert n%w == 1 and max(p,q) <= w*w and r > n//2
    assert n-r <= min(w,2*p-1,2*q-1)
    assert block_fails(n,p,w,w) and block_fails(n,q,w,w)
    assert q//p >= 3 and q < p*p and p >= 3
    assert (n-1)%2 == (n-1)%w == 0
    print(f'Optimal four-prime certificate at N={n} with fourth prime N-94.')


if __name__ == '__main__':
    exact_controls()
    square_root_controls()
    four_primes()
