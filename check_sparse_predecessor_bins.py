"""Check retained prime powers by factorial floors and, at small weights, exact coefficients."""

from collections import Counter
from functools import lru_cache
from itertools import combinations_with_replacement, product
from math import comb, gcd

from check_three_factors import powers, valuation


factor = lru_cache(None)(powers)


def check(base, row):
    N = sum(row)
    M = N-1
    bins = [1, 1, 1]
    J, G = 1, 1
    smooth = 1
    for prime, power in factor(M).items():
        value = valuation(N, row, prime)
        retained = power
        while value and retained > 1:
            value -= 1
            retained //= prime
        G *= power//retained
        if base % prime == 0:
            smooth *= power
            continue
        if retained > 1:
            remainders = [v % retained for v in row]
            assert sum(remainders) == 1
            bins[remainders.index(1)] *= retained
            J *= retained
    number = sum(v > 1 for v in bins)
    if number <= 1:
        assert sum(a % J == 0 for a in row) >= 2
    elif number == 2:
        assert J**5 < M**2
    else:
        assert J**35 < M**12
    if smooth == 1 or (len(factor(base)) == 1 and smooth**5 <= M**3):
        assert G > 1
    if N <= 1000:
        a,b,c = row
        K = comb(N, a)*comb(b+c, b)
        assert gcd(M,K) == G
    return number, N <= 1000


def sparse_rows():
    counts = Counter()
    for base in (2, 3, 4, 5, 6, 10, 15):
        atoms = [base**i for i in range(5)]
        sums = {n: sorted({sum(t) for t in combinations_with_replacement(atoms, n)})
                for n in range(1, 5)}
        # Sort lengths, then suppress only permutations of equal-length parts.
        for lengths in combinations_with_replacement(range(1,5), 3):
            if sum(lengths) > 6:
                continue
            for row in product(*(sums[t] for t in lengths)):
                if any(lengths[i]==lengths[i+1] and row[i]>row[i+1] for i in (0,1)):
                    continue
                occupied, exact = check(base, row)
                counts["rows"] += 1
                counts["exact products"] += exact
                counts[f"{occupied} occupied"] += 1
    return counts


def bin_shapes():
    # All three lowest-layer groups really can occur. Higher carries remove
    # the displayed primes from the retained product in this example.
    row = [178**i*(1+178**d) for i,d in ((3,15),(2,21),(3,35))]
    N = sum(row)
    for j, (prime,value) in enumerate(((29,29),(11,42),(7,52))):
        assert [a%prime for a in row] == [int(i==j) for i in range(3)]
        assert valuation(N,row,prime) == value
    checks = 0
    for u,v,w in combinations_with_replacement(range(3,40,2),3):
        if gcd(u,v)*gcd(v,w)*gcd(u,w) != 1:
            continue
        assert 5*(v+w-2) <= 2*v*w
        assert 35*(u+v+w-3) <= 12*v*w
        assert gcd(u*v,u*w) == u
        checks += 1
    return checks


if __name__ == "__main__":
    print(sparse_rows())
    print("coprime index triples", bin_shapes())
