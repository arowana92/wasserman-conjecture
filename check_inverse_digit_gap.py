"""Check the unique inverse pair forced by a second-interval obstruction."""

from collections import Counter
from itertools import combinations
from math import comb, gcd, prod

from check_endpoint_pair_bound import prime_powers


def adds_freely(x, y, p):
    while x and y:
        x, a = divmod(x, p)
        y, b = divmod(y, p)
        if a+b >= p:
            return False
    return True


def has_digit_gap(n, p, length):
    last = None
    i = 0
    while n:
        n, digit = divmod(n, p)
        if digit:
            if last is not None and i-last > length:
                return True
            last = i
        i += 1
    return False


def gap_controls():
    counts = Counter()
    for p in (2,3,5,7):
        for a in (1,2,3):
            power = p**a
            for q in range(2,1500):
                if q%p == 0:
                    continue
                u = pow(q,-1,power)
                v = (power*q+1-q*u)//power
                assert power*v + q*u == power*q+1 and 0 < v < q
                gap = has_digit_gap(q,p,a)
                if gap:
                    assert not adds_freely(v,q-v,p)
                if q <= 120:
                    assert adds_freely(v,q-v,p) == (comb(q,v)%p != 0)
                    counts['exact inverse binomials'] += 1
                counts['numerical digit gaps'] += gap
                counts['no gap, carry-free inverse'] += not gap and adds_freely(v,q-v,p)
    return counts


def inverse_pairs(limit):
    atoms = []
    for power in range(2,limit+1):
        factors = tuple(prime_powers(power))
        if len(factors) == 1:
            atoms.append(factors[0])
    counts = Counter()
    surviving_binary = []
    for (p,P),(q,Q) in combinations(atoms,2):
        if p == q:
            continue
        u, v = pow(Q,-1,P), pow(P,-1,Q)
        x, y = P*v, Q*u
        assert x+y == P*Q+1 and gcd(x,y) == gcd(u,v)
        counts['all pairs'] += 1
        bases_free = adds_freely(x,y,p) and adds_freely(x,y,q)
        admissible = prod(w for s,w in prime_powers(gcd(u,v)) if adds_freely(x,y,s))
        allowed = admissible >= 3 and bases_free
        binary = allowed and (x & y) == 0
        counts['passes inverse and scale tests'] += allowed
        counts['passes with binary marker'] += binary
        if binary:
            surviving_binary.append((P,Q))
        if P*Q <= 1200:
            inverse_rows = [(k*P,P*Q+1-k*P) for k in range(1,Q+1)
                            if (P*Q+1-k*P)%Q == 0]
            assert inverse_rows == [(x,y)]
            value = comb(x+y,x)
            assert bases_free == (gcd(p*q,value) == 1)
            decomp = tuple(prime_powers(gcd(x,y)))
            assert admissible == prod(w for s,w in decomp if value%s != 0)
            assert ((x & y) == 0) == (value%2 != 0)
            counts['exact binomial controls'] += 1
        assert not (P <= 4 and allowed)
        assert not (P == 5 and binary)
    return counts, surviving_binary


if __name__ == '__main__':
    print(dict(gap_controls()))
    counts, pairs = inverse_pairs(1000)
    print(dict(counts))
    print('Binary pairs left by these necessary tests:', pairs)
