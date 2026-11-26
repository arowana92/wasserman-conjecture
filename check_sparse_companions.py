"""Independently audit the remaining rows by single-coordinate CRT and factorial floors."""

from argparse import ArgumentParser
from collections import Counter
from math import prod

from flint import fmpz

from check_three_factors import valuation


def factor(n):
    entries = [(int(p), int(p)**e) for p, e in fmpz(n).factor()]
    if prod(w for p, w in entries) != n or any(not fmpz(p).is_prime() for p, w in entries):
        raise ArithmeticError(n)
    return entries


def eligible_entries(weight, target, T, atoms):
    # alpha = 0 or target modulo every w. Build the complete set one w at a time.
    entries, modulus = [0], 1
    for _, w in atoms:
        inv = pow(modulus, -1, w)
        entries = [a + modulus*(((b-a)*inv) % w) for a in entries for b in (0,target)]
        modulus *= w
    assert modulus == T and T > target
    scale = weight//target
    valid = sorted(x for x in entries if 0 < x < target)
    assert all((scale*x)*(scale*x-1) % T == 0 for x in valid)
    return valid


def inspect(p, P, n, S, counts):
    V = n//P
    assert V*P == n and V % p and factor(P) == [(p,P)] and (n-1) % S == 0
    T = (n-1)//S
    shifted = factor(n-1)
    needed = [(q,w) for q,w in shifted if S % q != 0]
    assert prod(w for q,w in needed) == T
    if P == p:
        assert V == 3*p+2
        first = (p*p,p*(p+1),p*(p+1))
        assert valuation(n,first,p) == 0
        assert S == prod(w for q,w in shifted if valuation(n,first,q) == 0)
    entries = eligible_entries(n, P, T, needed)
    counts['audited systems'] += 1
    counts['single entries'] += len(entries)
    seen = set(entries)
    for i, a in enumerate(entries):
        for b in entries[i:]:
            c = P-a-b
            if c < b:
                break
            if c not in seen:
                continue
            counts['positive triples of entries'] += 1
            row = (V*a,V*b,V*c)
            if any(valuation(n,row,q) for q,_ in needed):
                continue
            counts['passing all predecessor primes'] += 1
            if all(valuation(n,row,q) == 0 for q,_ in factor(V)):
                raise ArithmeticError((p,P,n,S,row))


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('logs', nargs='+')
    args = parser.parse_args()
    counts = Counter()
    for path in args.logs:
        with open(path) as stream:
            systems = [tuple(map(int,line.split()[1:])) for line in stream if line.startswith('CRT ')]
        assert len(systems) == len(set(systems))
        for data in systems:
            inspect(*data, counts)
    print(dict(counts))


if __name__ == '__main__':
    main()
