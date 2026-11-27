"""Independent all-exponent check using rational bounds for the logarithm series."""

from argparse import ArgumentParser
from fractions import Fraction
from functools import lru_cache
from math import isqrt


S = 10**70
F = S * 10**8


def ceil_div(n, d):
    return -(-n//d)


@lru_cache(maxsize=260)
def log_interval(n, d):
    """Outward fixed-point enclosure of ln(n/d) for 1 <= n/d <= 2."""
    a, b = n-d, n+d
    assert 0 <= 3*a <= b
    lower, upper = F*a//b, ceil_div(F*a,b)
    lo = hi = 0
    j = 1
    while upper > 2:
        lo += lower//j
        hi += ceil_div(upper,j)
        lower = lower*a*a//(b*b)
        upper = ceil_div(upper*a*a,b*b)
        j += 2
    # The omitted positive series before multiplying by 2 is <= 2*u^j.
    return 2*lo, 2*hi+4*upper


def binary_interval(x):
    a,b = x.numerator,x.denominator
    e = 0
    while a < b:
        a *= 2
        e -= 1
    while a >= 2*b:
        b *= 2
        e += 1
    k = 256*a//b
    v,w = log_interval(k,256)
    v1,w1 = log_interval(256*a,k*b)
    lo2,hi2 = log_interval(2,1)
    return e*S+(v+v1)*S//hi2, e*S+ceil_div((w+w1)*S,lo2)


def nearest(word):
    return abs((word+S//2)%S-S//2)


def prime_list(limit):
    # Remove odd composite indices 2j+1 by marking only odd multiples.
    flags = bytearray([1])*((limit+1)//2)
    flags[0] = 0
    for j in range(1,(isqrt(limit)+1)//2):
        if flags[j]:
            p = 2*j+1
            first = (p*p-1)//2
            flags[first::p] = bytes((len(flags)-1-first)//p+1)
    return [2*i+1 for i,flag in enumerate(flags) if flag]


def run(limit):
    fractions = {4587:[Fraction(256,417)], 7383:[Fraction(19597,35328), Fraction(486233,876544),
                Fraction(836139,1507328), Fraction(46678369,84148224)], 9453:[Fraction(256,411)]}
    logs = []
    for D,limits in fractions.items():
        for f in limits:
            for c in (8*D*f/5,160*f/D):
                logs.append(binary_interval(c))
    primes = prime_list(limit)
    largest = 0
    for p in primes:
        lo,hi = binary_interval(Fraction(p))
        unresolved = list(logs)
        denominator, previous = 0, 1
        top, bottom = lo, S
        while bottom and unresolved:
            q,r = divmod(top,bottom)
            previous,denominator = denominator,q*denominator+previous
            top,bottom = bottom,r
            if denominator <= 10**16:
                continue
            assert denominator < 10**32
            loss = 6*10**14*(nearest(denominator*lo)+denominator*(hi-lo))
            unresolved = [(a,b) for a,b in unresolved
                          if 256*(nearest(denominator*a)-denominator*(b-a)-loss) <= S]
        assert not unresolved, p
        largest = max(largest, denominator)
    print(f'{len(primes)} odd bases, {len(primes)*len(logs)} certified offsets, '
          f'maximum final denominator {largest}.')


if __name__ == '__main__':
    cli = ArgumentParser(description=__doc__)
    cli.add_argument('bound',type=int,nargs='?',default=1000000)
    args = cli.parse_args()
    assert 3 <= args.bound <= 1000000
    run(args.bound)
