"""Reconstruct all endpoint assignments of a prescribed product, independent of weight."""

from argparse import ArgumentParser
from collections import Counter
from functools import lru_cache
from itertools import permutations
from math import gcd, prod

from check_cross_parts import divisors
from check_three_factors import valuation
from search_even_endpoint import prime_powers


@lru_cache(maxsize=32768)
def factor_pairs(kernel):
    return tuple((v, kernel//v) for v in divisors(kernel))


def numerators(n, slots, lower=1):
    if slots == 1:
        if n >= lower:
            yield (n,)
        return
    for v in range(lower, n+1):
        if v**slots > n:
            break
        if n % v == 0:
            for tail in numerators(n//v, slots-1, v):
                yield (v,) + tail


def fraction_data(coeffs, bins):
    product = prod(bins)
    numerator = sum(c*(product//w) for c,w in zip(coeffs,bins))
    return product, numerator


def assignments(first, second):
    total = sum(first)*sum(second)
    found = set()

    def extend(s, t):
        i, j = len(s), len(t)
        S, X = fraction_data(first, s)
        T, Y = fraction_data(second, t)
        K = S*T
        deficit = K-X*Y
        if deficit <= 0:
            return
        if i+j < 4:
            # The next bin of at least one side lies below its displayed cap.
            left = total-sum(first[:i])*sum(second[:j])
            if i < 3:
                cap = (left*K-1)//(deficit*t[0])
                for w in range(s[-1], cap+1):
                    if gcd(w, K) == 1:
                        extend(s+(w,), t)
            if j < 3:
                cap = (left*K-1)//(deficit*s[0])
                for w in range(t[-1], cap+1):
                    if gcd(w, K) == 1:
                        extend(s, t+(w,))
            return

        # The two unfilled bins satisfy deficit*v*w - c*w - e*v = constant.
        missing = [(0, f*S*Y) for f in first[i:]] + [(1, g*T*X) for g in second[j:]]
        (side_v, c), (side_w, e) = missing
        kernel = c*e-deficit
        if side_v != side_w:
            kernel = first[-1]*second[-1]*K*K-deficit
        if kernel <= 0:
            return
        for x, y in factor_pairs(kernel):
            v, remainder = divmod(x+c, deficit)
            w, other_remainder = divmod(y+e, deficit)
            if remainder or other_remainder:
                continue
            additions = ((side_v, v), (side_w, w))
            s1 = s+tuple(b for side,b in additions if side == 0)
            t1 = t+tuple(b for side,b in additions if side == 1)
            if tuple(sorted(s1)) != s1 or tuple(sorted(t1)) != t1:
                continue
            ss, V = fraction_data(first, s1)
            tt, U = fraction_data(second, t1)
            assert U*V == ss*tt+1
            if gcd(ss, tt) != 1 or gcd(U, V) != 1:
                continue
            rows = (tuple(sorted(zip(first,s1))), tuple(sorted(zip(second,t1))))
            found.add(tuple(sorted(rows)))

    for s in range(1, total):
        for t in range(1, (total-1)//s+1):
            if gcd(s, t) == 1:
                extend((s,), (t,))
    return found


def inspect(value):
    systems = set()
    for unordered in numerators(value, 6):
        for order in set(permutations(unordered)):
            first, second = order[:3], order[3:]
            if sorted(first) <= sorted(second):
                systems.update(assignments(first, second))
    counts = Counter()
    largest = 0
    for system in systems:
        normalized = []
        for entries in system:
            modulus = prod(w for _, w in entries)
            normalized.append([c*(modulus//w) for c,w in entries])
        a, b = normalized
        n = sum(a)*sum(b)
        largest = max(largest,n)
        row = [x*sum(b) for x in a]
        other = [y*sum(a) for y in b]
        assert prod(row)*prod(other) == value*n**3*(n-1)**2
        counts['assignments'] += 1
        if all(r[0] | r[1] | r[2] != n for r in (row, other)):
            continue
        counts['after binary'] += 1
        for p in sorted(set(prime_powers(n)) | set(prime_powers(n-1))):
            if valuation(n,row,p) and valuation(n,other,p):
                counts['common prime found'] += 1
                break
        else:
            raise ArithmeticError(('requires independent counterexample check', n, row, other))
    print(value, dict(counts), 'largest weight', largest, flush=True)
    return systems


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('limit', type=int, nargs='?', default=32)
    args = parser.parse_args()
    for value in range(2, args.limit+1, 2):
        inspect(value)
