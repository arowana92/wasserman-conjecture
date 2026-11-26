"""Check the two predecessor-power bound and its leading-power consequences."""

from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations
from math import gcd, prod

from check_endpoint_pair_bound import prime_powers, rows
from check_scaled_row_blocks import certify_prime


def check_pair_bounds(bound):
    counts = Counter()
    for n in range(3, bound+1):
        table = list(rows(n))
        for (_, p), (_, q) in combinations(prime_powers(n-1), 2):
            for parts, first in table:
                smallest = sum(parts[:2])
                for _, second in table:
                    shared = gcd(first, second)
                    divisor = gcd(shared, n*p*q)
                    assert divisor*smallest > p*q
                    counts['ordered row pairs and prime-power choices'] += 1
                    counts['pairs requiring a common divisor by size'] += smallest <= p*q
    return counts


def subpowers(prime, full):
    power = prime
    while power <= full:
        yield power
        power *= prime


def check_certificates(bound):
    counts = Counter()
    for n in range(3, bound+1):
        weight = {s for s, _ in prime_powers(n)}
        atoms = tuple(prime_powers(n-1))
        values = [value for _, value in rows(n)]
        seen = set()
        for (p, full_p), (q, full_q) in combinations(atoms, 2):
            for s in weight | {a for a, _ in atoms}:
                leading = s
                while s*leading <= n:
                    leading *= s
                if 2*leading <= n:
                    continue
                off = n-leading
                for power_p in subpowers(p, full_p):
                    for power_q in subpowers(q, full_q):
                        product = power_p*power_q
                        if off > 2*product:
                            continue
                        if off > product and gcd(n, product+1) > 2 and len(tuple(prime_powers(product+1))) != 1:
                            continue
                        counts['second-interval parameter choices' if off > product else 'first-interval parameter choices'] += 1
                        prime_list = tuple(sorted(weight | {p, q, s}))
                        if prime_list in seen:
                            continue
                        seen.add(prime_list)
                        counts['distinct weight and prime-list choices'] += 1
                        counts['choices omitting a predecessor prime'] += any(a not in prime_list for a, _ in atoms)
                        masks = {sum(1 << i for i, a in enumerate(prime_list) if value%a == 0) for value in values}
                        assert all(x & y for x in masks for y in masks), (n, prime_list)
    return counts


def example():
    p, q = 2305843009213693951, 147573952589676412931
    factors = (2, 3, 3, 307, 2857, 6529, 43691)
    for prime in {p, q, *factors}:
        certify_prime(prime)
    r = prod(factors)
    n = 1+r*p*q
    assert r == 4503599627370498 and n == 1532495540865889538289620652285986858360875510602399739
    off = n-2**180
    assert n.bit_length() == 181 and p*q < off < 2*p*q
    assert off**3 > n*n and bin(off).count('1') == 116 and gcd(n, p*q+1) == 1
    print('181-bit example verified.')


def subset_sums(bits):
    values = [0]
    for bit in bits:
        values += [x+bit for x in values]
    return values


def rotation_windows():
    threshold = Fraction(385, 383)
    for p, q, first_cover in ((2,3,253), (3,2,401)):
        for length in (first_cover-1, first_cover):
            positions = []
            point = Fraction(1)
            for _ in range(length):
                positions.append(point)
                point *= q
                while point >= p:
                    point /= p
            positions.sort()
            ends = positions[1:]+[p*positions[0]]
            gap = max(y/x for x,y in zip(positions, ends))
            if length < first_cover:
                assert gap == Fraction(2**233,3**147) and gap > threshold
            else:
                assert gap == Fraction(2**317,3**200) and gap < threshold
    print('Rotation windows verified with exact rational quotients.')


def four_prime_certificates():
    n = 1208925867235455872460352
    weight_prime = n//64
    near_prime = n-35
    p, q = 134217689, 268947953
    for s in (p, q, weight_prime, near_prime, 283, 487):
        certify_prime(s)
    assert n == 64*weight_prime and n-1 == 3**5*283*487*p*q
    offset = n-2**80
    assert p*q < offset < 2*p*q and gcd(n, p*q+1) == 2
    assert 35 < 64 < p < q < weight_prime and 2*near_prime > n
    bits = [1 << i for i in range(offset.bit_length()) if (offset >> i) & 1]
    assert len(bits) == 29
    left, right = subset_sums(bits[:14]), subset_sums(bits[14:])
    for prime, expected in ((p, (4, 0)), (q, (4, 1))):
        residues = defaultdict(list)
        for x in left:
            residues[x%prime].append(x)
        candidates = [[x+y for y in right for x in residues[(a-y)%prime] if x+y]
                      for a in (0, 1)]
        assert tuple(map(len, candidates)) == expected
        assert all(x & y for x in candidates[0] for y in candidates[0]+candidates[1])
    print('Two optimal four-prime certificates verified at', n)


if __name__ == '__main__':
    print(dict(check_pair_bounds(48)))
    print(dict(check_certificates(240)))
    example()
    rotation_windows()
    four_prime_certificates()
