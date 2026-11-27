"""Reconstruct cuts from divisors of the part itself and verify the native sieve."""

from argparse import ArgumentParser
from collections import Counter
from functools import lru_cache
from itertools import combinations, combinations_with_replacement, product
from math import comb, gcd, lcm, prod
from pathlib import Path
import random
import re
import subprocess

from flint import fmpz


@lru_cache(maxsize=None)
def factor(n):
    entries = [(int(p), int(e)) for p, e in fmpz(n).factor()]
    assert prod(p**e for p, e in entries) == n and all(fmpz(p).is_prime() for p, _ in entries)
    return entries


def divisors(fs):
    return [prod(values) for values in product(*([p**i for i in range(e+1)] for p, e in fs))]


def value(n, row, p):
    total = 0
    layer = p
    while layer <= n:
        total += n//layer - sum(z//layer for z in row)
        layer *= p
    return total


def assignments(n, retained, stats):
    nf = factor(n)
    pred = sorted(factor(n-1), key=lambda pair: pair[0]**pair[1], reverse=True)
    weights = divisors(nf) if retained else divisors([(p**e, 1) for p, e in nf])
    sizes = [n//u for u in weights if u >= (2 if retained else 3) and n//u >= 3]
    for exponents in product(*(range(e+1) if retained else (0, e) for _, e in pred)):
        atoms = [(p**b, p, e-b) for (p, e), b in zip(pred, exponents) if b]
        if len(atoms) < 3:
            continue
        P = prod(w for w, _, _ in atoms)
        for v in sizes:
            if not P//atoms[-1][0] < v <= P or 2*max(w for w, _, _ in atoms)+1 > v:
                continue
            u = n//v
            stats[0] += 1
            cuts = {}
            indexes = frozenset(range(len(atoms)))
            for r in range(1, len(atoms)-1):
                for subset in combinations(indexes, r):
                    b = prod(atoms[i][0] for i in subset)
                    multiple = u*(P//b)
                    a = multiple*pow(multiple, -1, b)
                    if a >= n:
                        continue
                    cuts[frozenset(subset)] = a
            # The native cut statistic includes masks leaving just one atom,
            # although they cannot participate in a positive three-way split.
            for i in indexes:
                b = P//atoms[i][0]
                multiple = u*atoms[i][0]
                a = multiple*pow(multiple, -1, b)
                if a < n:
                    cuts[indexes-{i}] = a
            for group, a in cuts.items():
                if 3*a > n or any(value(n, (a, n-a), p) > loss for _, p, loss in atoms):
                    continue
                stats[1] += 1
                left = indexes-group
                for second, b in cuts.items():
                    if not second < left or b < a:
                        continue
                    c = cuts.get(left-second, 0)
                    if c < b or a+b+c != n:
                        continue
                    row = (a, b, c)
                    if all(value(n, row, p) <= loss for _, p, loss in atoms):
                        yield n, u, P, *row


def compare(binary, weights, retained):
    lines = []
    stats = [0, 0]
    expected = Counter()
    for n in weights:
        data = [n]
        for m in (n, n-1):
            fs = factor(m)
            data += [len(fs), *(v for p, e in fs for v in (p, p**e))]
        lines.append(" ".join(map(str, data)))
        expected.update(assignments(n, retained, stats))
    command = [str(Path(binary).resolve())] + (["--retained"] if retained else [])
    result = subprocess.run(command, input="\n".join(lines)+"\n", capture_output=True,
                            text=True, check=True)
    summary = list(map(int, re.findall(r"\d+", result.stderr)))
    found = Counter()
    for line in result.stdout.splitlines():
        n, u, P, a, b, c = map(int, line.split())
        found[(n, u, P, *sorted((a, b, c)))] += 1
    assert expected == found and summary == [len(weights), *stats, sum(found.values())]
    if 222788354 in weights:
        assert (222788354, 283, 887603, 45491401, 86617810, 90679143) in found
        assert sum(k[0] == 8284703 for k in found) >= 1
    print("Retained" if retained else "Full", summary)


def exact_controls():
    tests = 0
    for n in range(3, 161):
        for a in range(1, n//3+1):
            for b in range(a, (n-a)//2+1):
                row = (a, b, n-a-b)
                k = comb(n, a)*comb(n-a, b)
                g = gcd(*row)
                assert gcd(n-1, k) >= g
                for p, e in factor(n-1):
                    valuation = value(n, row, p)
                    assert k % p**valuation == 0 and k % p**(valuation+1) != 0
                    w = p**max(0, e-valuation)
                    if w > 1:
                        assert sorted(x%w for x in row) == [0, 0, 1]
                tests += 1
    return tests


def subdivisor_controls():
    checked = 0
    for length in range(3, 7):
        for row in combinations_with_replacement(range(1, 9), length):
            n = sum(row)
            remaining, coefficient = n, 1
            for a in row:
                coefficient *= comb(remaining, a)
                remaining -= a
            U = gcd(*row)
            for D in divisors(factor(n-1)):
                if D == 1:
                    continue
                h = len(factor(D))
                if h < length:
                    assert gcd(D, coefficient)*(n//U-h) >= (length-h)*D
                    checked += 1
    return checked


def transfer_controls():
    checked = nontrivial = at_base = 0
    for d in range(2, 19):
        period = lcm(*(p**(e-1)*(p-1) for p, e in factor(lcm(*range(1, d+2)))))
        for q in (2, 3, 5, 7, 19):
            Q = q**period
            n = 1+d*Q
            for a in range(1, d//2+1):
                for b in range(a, d-a+1):
                    c = d-a-b
                    if value(d, (a, b, c), q):
                        continue
                    row = (a*Q, b*Q, c*Q+1)
                    U = gcd(*row)
                    small = (a, b, c+1)
                    assert gcd(*small) % U == 0
                    small_J = gcd(d, comb(d+1, a)*comb(d+1-a, b))
                    large_J = prod(p**min(e, value(n, row, p)) for p, e in factor(d))
                    assert large_J % small_J == 0 and large_J >= U
                    checked += 1
                    nontrivial += U > 1
                    at_base += d % q == 0
    n = 222788354
    row = (45491401, 86617810, 90679143)
    endpoints = [prod(p**min(e, value(n, row, p)) for p, e in factor(m)) for m in (n, n-1)]
    assert gcd(*row) == 283 and endpoints == [787238, 251]
    return checked, nontrivial, at_base


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("binary")
    args = parser.parse_args()
    print("Exact multinomial controls:", exact_controls(), flush=True)
    print("Partial predecessor divisor bounds:", subdivisor_controls())
    print("Transferred rows, nontrivial part gcds, base dividing cofactor:", transfer_controls())
    weights = set(m+1 for m in divisors(factor(lcm(*range(1, 24)))) if m > 1)
    rng = random.Random(222788354)
    primes = [int(p) for p, _ in factor(lcm(*range(1, 62)))]
    for _ in range(90):
        m = 1
        for p in rng.sample(primes, 6):
            w = p
            while p*w <= 61 and rng.randrange(2):
                w *= p
            m *= w
        weights.add(m+1)
    weights.update((8284703, 222788354, 1+100003*100019*100043*100049))
    for retained in (False, True):
        compare(args.binary, sorted(weights), retained)
