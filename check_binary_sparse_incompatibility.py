"""Audit joint base avoidance by exact coefficients and by digit-unit assignments."""

from collections import Counter, defaultdict
from itertools import combinations, product
from math import comb, gcd, isqrt


def units(n, base):
    values = []
    exponent = 0
    while n:
        n, digit = divmod(n, base)
        values.extend([exponent]*digit)
        exponent += 1
    return values


def v2(n):
    return (n & -n).bit_length()-1


def coefficient(row):
    a, b, c = row
    return comb(a+b+c, a)*comb(b+c, b)


def primes(limit):
    flags = bytearray([1])*(limit+1)
    flags[:2] = bytes(2)
    for p in range(2, isqrt(limit)+1):
        if flags[p]:
            flags[p*p::p] = bytes((limit-p*p)//p+1)
    return [p for p, flag in enumerate(flags) if flag]


def binary_bounds():
    counts = Counter()
    odd = primes(53)[1:]
    for n in range(3, 321):
        bits = units(n, 2)
        for a in range(1, n//3+1):
            for b in range(a, (n-a)//2+1):
                c = n-a-b
                if a&b or (a|b)&c:
                    continue
                k = coefficient((a,b,c))
                assert k % 2
                counts['binary-free exact coefficients'] += 1
                for q in odd:
                    if k % q == 0:
                        continue
                    eta = n%2
                    required = eta + 4*(3-eta) - 2*len(set(bits) & {1, v2(q+1)})
                    second = sum(2**min(e, v2(q-1)) for e in bits[:3])
                    assert len(units(n,q)) >= max(required, second)
                    counts['joint avoidances'] += 1
    for row, q, cost in [((12,162,768),3,8), ((80,130,260),5,10)]:
        assert sum(row) % q == 0 and len(units(sum(row),q)) == cost
        assert gcd(coefficient(row), 2*q) == 1
    print(dict(counts))


def check_certificate():
    n = 396
    members = (2,3,11,389)
    counts = Counter()
    for a in range(1,n//3+1):
        for b in range(a,(n-a)//2+1):
            k = coefficient((a,b,n-a-b))
            avoided = [p for p in members if k % p]
            assert len(avoided) <= 1
            counts[tuple(avoided)] += 1
    assert all(counts[(p,)] > 0 for p in members)
    print(f'{sum(counts.values())} exact coefficients at 396, avoidance counts {dict(counts)}.')


def even_order(a, p):
    order = p-1
    while order % 2 == 0 and pow(a,order//2,p) == 1:
        order //= 2
    return order % 2 == 0


def check_weight_pairs(limit):
    sparse = defaultdict(list)
    for p in primes(limit//2)[1:]:
        for n in range(2*p,limit+1,2*p):
            positions = units(n,p)
            if len(positions) <= 6:
                sparse[n].append((p,positions))
    counts = Counter()
    for n, bases in sparse.items():
        for (p, pos), (q, other) in combinations(bases, 2):
            obstruction = (len(pos) < 6 or len(other) < 6 or n%4 == 0
                           or not (even_order(q,p) and even_order(p,q)))
            if not obstruction:
                continue
            counts['incompatible weight-prime pairs'] += 1
            pieces = [p**e for e in pos]
            # Assign the (possibly repeated) units individually, with all three slots used.
            for choice in product(range(3), repeat=len(pieces)):
                if set(choice) != {0,1,2}:
                    continue
                row = [sum(x for x, slot in zip(pieces,choice) if slot == j) for j in range(3)]
                assert sum(len(units(x,q)) for x in row) > len(other)
                counts['ordered digit allocations'] += 1
    print(dict(counts))


if __name__ == '__main__':
    binary_bounds()
    check_certificate()
    check_weight_pairs(200_000)
