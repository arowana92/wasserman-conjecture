"""Recover the rational binary limits and certify all exponents at bounded prime bases."""

from argparse import ArgumentParser
from fractions import Fraction
from math import isqrt

from check_endpoint_pair_bound import prime_powers
from check_seven_quadratic_suffixes import RESIDUAL, upper_graph

BITS = 200
SCALE = 1 << BITS
MAX_EXPONENT = 10**18


def limits(graph):
    live = set(range(len(graph)))
    removed = 0
    while True:
        new = {i for i in live if any(target in live for _,target in graph[i])}
        if new == live:
            break
        live = new
        removed += 1
    todo = [(0,1,1,{})]
    values = []
    while todo:
        i,word,h,seen = todo.pop()
        if i in seen:
            length, prefix = seen[i]
            period = h-length
            cycle = word-(prefix<<period)
            for node,(at,_) in seen.items():
                if at >= length:
                    assert sum(dst in live for _,dst in graph[node]) == 1
            values.append(Fraction(prefix*(2**period-1)+cycle, 2**length*(2**period-1)))
            continue
        previous = {**seen, i:(h,word)}
        for bit,target in graph[i]:
            if target in live:
                todo.append((target,2*word+bit,h+1,previous))
    return removed, sorted(values)


def constants():
    rays = []
    for D,a,b,z,w in RESIDUAL:
        _,graph = upper_graph(D,z,w)
        dead, values = limits(graph)
        assert dead == (6 if D == 7383 else 0)
        print(f'D={D}: maximum trimmed tail {dead}, binary limits ' + ', '.join(map(str,values)))
        for v in values:
            c = Fraction(8*D,5)*v
            rays.extend((c, 100*c/(D*D)))
    assert len(rays) == 12
    assert all(max(c.numerator,c.denominator)<10**11 and Fraction(1,1000)<c<10000 for c in rays)
    for c in rays:
        support = {p for v in (c.numerator,c.denominator) for p,_ in prime_powers(v) if p != 2}
        assert len(support) >= 2
    assert 26*2**38*38 + 16 < 3*10**14
    return rays


def normalize(x):
    n, d = x.numerator, x.denominator
    e = n.bit_length()-d.bit_length()
    if e >= 0:
        d <<= e
    else:
        n <<= -e
    if n < d:
        n *= 2
        e -= 1
    assert d <= n < 2*d
    return n,d,e


def log_floor(x):
    n,d,e = normalize(x)
    precision = 3*BITS+32
    while True:
        F = 1 << precision
        lo, hi = (F*n)//d, (F*n+d-1)//d
        word = 0
        for _ in range(BITS):
            lo = lo*lo >> precision
            hi = (hi*hi+F-1) >> precision
            if lo < 2*F <= hi:
                break
            bit = int(lo >= 2*F)
            if bit:
                lo //= 2
                hi = (hi+1)//2
            word = 2*word+bit
        else:
            return e*SCALE+word
        precision *= 2


def primes(bound):
    flags = bytearray([1])*(bound+1)
    for n in range(2, isqrt(bound)+1):
        if flags[n]:
            flags[n*n::n] = bytes((bound-n*n)//n+1)
    return (p for p in range(3,bound+1,2) if flags[p])


def distance(word):
    residue = word%SCALE
    return min(residue,SCALE-residue)


def denominators(A):
    a, b = A, SCALE
    old, current = 1, 0
    while b:
        q, remainder = divmod(a,b)
        old,current = current,q*current+old
        if current >= 10**32:
            break
        if current > MAX_EXPONENT:
            yield current
        a,b = b,remainder


def run(bound):
    rays = constants()
    offsets = [log_floor(c) for c in rays]
    total = maximum = 0
    for p in primes(bound):
        A = log_floor(Fraction(p*p))
        pending = list(offsets)
        for v in denominators(A):
            loss = v + MAX_EXPONENT*(distance(v*A)+v)
            pending = [B for B in pending if 1000*(distance(v*B)-loss) <= SCALE]
            if not pending:
                maximum = max(maximum,v)
                break
        assert not pending, (p,pending)
        total += 1
    print(f'{total} odd prime bases through {bound}; {total*len(rays)} logarithmic offsets '
          f'separated; largest final denominator {maximum}.')


if __name__ == '__main__':
    cli = ArgumentParser(description=__doc__)
    cli.add_argument('bound',type=int,nargs='?',default=1000000)
    args = cli.parse_args()
    if not 3 <= args.bound <= 1000000:
        cli.error('bound must lie in [3,1000000]')
    run(args.bound)
