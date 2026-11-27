"""Exact certificates for joint_prime_counts.md. Uses only the standard library."""

from fractions import Fraction as F
from functools import lru_cache
from math import comb, isqrt, prod

from check_prime_residue_counts import log_interval, log_up


PRIMES = (2, 3, 5, 7, 11)
TARGETS = {3: F(7593, 10000), 4: F(7265, 10000), 5: F(7106, 10000)}

# A mask assigns its set bits to the smaller product. Each line lists the
# mask, a, b, the complement tilts and the row tilts, in units of 1/100000
# except for the mask. Blocks correspond to 3, 4 and 5 successive primes.
CERTIFICATES = """
0 1 74999 56378 64968 74908 56378 64968 74908
1 1 74999 99998 47140 60127 21853 47139 60126
2 1 74999 38287 99998 61189 38286 31764 61188
3 59553 15230 48843 58184 86457 40242 50259 48535
4 1 74999 40256 50246 99999 40256 50246 46302
5 51026 23915 52076 72827 71683 38239 36906 61232
6 45647 29477 60585 63102 73282 28118 47107 60156
7 74999 1 56378 64968 74908 56378 64968 74908

0 1 74999 65142 72465 80631 85067 65142 72465 80631 85067
1 1 74999 99998 62582 73034 78938 27243 62581 73034 78938
2 1 74999 54412 99999 73556 79365 54412 37267 73556 79365
3 4115 70071 95279 96409 59994 67886 40251 50250 56855 65176
4 1 74999 55578 64268 99999 80022 55577 64267 51221 80022
5 1 74999 99999 45644 99999 66873 38284 45643 61191 66872
6 1 74999 36802 99999 99999 67848 36801 47138 60127 67847
7 74999 1 56379 64968 74909 99999 56378 64968 74908 55795
8 1 74999 56379 64968 74909 99999 56378 64968 74908 60257
9 1 74999 99999 47140 60127 99999 36799 47139 60126 67849
10 1 74999 38287 99999 61189 99999 38286 45637 61188 66878
11 74999 1 55578 64268 99999 80022 55577 64267 46470 80022
12 1 74999 40257 50247 99999 99999 40256 50246 56861 65170
13 71235 3624 56275 94929 74822 80393 54405 33861 73558 79368
14 66924 7822 86791 66091 75763 81151 25354 62573 73035 78943
15 74999 1 65142 72465 80631 85067 65142 72465 80631 85067

0 1 74999 70493 76907 83915 87663 91587 70493 76907 83915 87663 91587
1 1 74999 99998 70341 79033 83793 88857 29086 70340 79033 83793 88857
2 1 74999 63129 99999 79357 84052 89041 63129 39150 79357 84051 89041
3 1 74999 99999 99999 70558 76899 83864 40255 50246 70557 76898 83863
4 1 74999 63914 71432 99999 84450 89324 63913 71432 52884 84450 89324
5 1 74999 99999 60693 99999 77702 84455 38284 60692 61190 77702 84455
6 1 74999 52347 99999 99999 78175 84802 52346 47139 60126 78175 84802
7 23651 49805 84343 87977 91782 75852 82983 56364 64967 74920 64280 74211
8 1 74999 64452 71886 80197 99999 89516 64452 71885 80197 61692 89516
9 1 74999 99999 61522 72193 99999 84855 36798 61522 72192 67850 84855
10 1 74999 53252 99999 72745 99999 85187 53251 45641 72744 66874 85186
11 18241 55161 87429 90386 67128 95036 81558 55562 64264 56182 80035 74461
12 1 74999 54469 63293 99999 99999 85695 54468 63292 56861 65170 85694
13 10994 62481 91931 51973 95846 96859 79712 54396 43454 73560 79373 74918
14 6793 67192 40165 96002 97302 97963 78730 34662 62574 73035 78942 75747
15 74999 1 65142 72466 80631 85067 99999 65142 72465 80631 85067 68529
16 1 74999 65142 72466 80631 85067 99999 65142 72465 80631 85067 72168
17 1 74999 99999 62582 73034 78938 99999 34703 62581 73034 78938 75716
18 1 74999 54413 99999 73556 79365 99999 54412 43504 73556 79365 74881
19 11216 62240 91796 93760 64312 71569 97851 54454 63288 56848 65182 85705
20 1 74999 55578 64268 99999 80022 99999 55577 64267 54791 80022 73406
21 4580 69711 96326 49200 98128 69524 99053 53241 45626 72746 66888 85191
22 1117 73700 37707 99285 99520 68501 99758 36794 61520 72192 67853 84856
23 74999 1 64452 71886 80197 99999 89516 64452 71885 80197 57312 89516
24 1 74999 56379 64969 74909 99999 99999 56378 64968 74908 61877 72236
25 549 74360 99529 47569 60499 99821 99880 52344 47138 60127 78175 84803
26 1 74999 38287 99999 61189 99999 99999 38286 60684 61188 77706 84456
27 74999 1 63914 71432 99999 84450 89324 63913 71431 48166 84450 89324
28 1 74999 40257 50246 99999 99999 99999 40256 50246 70553 76897 83867
29 74999 1 63129 99999 79357 84052 89041 63129 34346 79357 84051 89041
30 74999 1 99998 70341 79033 83793 88857 24642 70340 79033 83793 88857
31 74999 1 70493 76907 83915 87663 91587 70493 76907 83915 87663 91587
"""


@lru_cache(maxsize=None)
def digit_cost(p, numerator):
    w = F(numerator, 100000)
    assert 0 < w <= 1
    return log_up(sum(w**d for d in range(p)), p), (p-1)*log_up(1/w, p)


def exponent(primes, mask, numbers):
    n = len(primes)
    assert len(numbers) == 2+2*n
    a, b = (F(x, 100000) for x in numbers[:2])
    c = 1-a-b
    assert min(a, b, c) > 0
    complement = [digit_cost(p, w) for p, w in zip(primes, numbers[2:2+n])]
    parts = [digit_cost(p, w) for p, w in zip(primes, numbers[2+n:])]
    small = [i for i in range(n) if mask & (1 << i)]
    large = [i for i in range(n) if i not in small]

    def delta(costs, indices):
        return max(0, 1+sum(costs[i][0]-1 for i in indices))

    baseline = a*delta(parts, small) + b*delta(parts, large)
    baseline += c*delta(complement, range(n))
    # Maximize the affine expression over every normalized digit sum in [0,1].
    for i in range(n):
        weight = a if i in small else b
        baseline += max(weight*parts[i][1]/3, c*complement[i][1])
    return (baseline+2*a/3)/(1+2*a/3)


def certify_assignments():
    blocks = CERTIFICATES.strip().split('\n\n')
    assert len(blocks) == len(TARGETS)
    for (n, target), block in zip(TARGETS.items(), blocks):
        lines = [list(map(int, line.split())) for line in block.splitlines()]
        assert [line[0] for line in lines] == list(range(1 << n))
        bounds = [exponent(PRIMES[:n], mask, numbers) for mask, *numbers in lines]
        upper = max(bounds)
        assert upper < target
        modulus = prod(PRIMES[:n])
        density = F(1 << n, modulus)
        assert sum(x*(x-1) % modulus == 0 for x in range(modulus)) == 1 << n
        print(f'Modulus {modulus}: {len(bounds)} assignments, density {density}, '
              f'largest upper bound {float(upper):.8f} (mask {bounds.index(upper)}).')


def threshold(t):
    assert F(2, 3) < t < 1
    return 2*t*t/(2*t-1) + (3*t-1)**2/(2*(3*t-2))


def certify_primorials():
    examples = {79: F(7, 10), 283: F(17, 25), 2113: F(67, 100)}
    total = F(0)
    # Round every reciprocal log down to the same rational grid.
    unit = 10**8
    for p in range(2, max(examples)+1):
        if any(p % j == 0 for j in range(2, isqrt(p)+1)):
            continue
        upper_log = log_interval(F(p))[1]
        total += F(int(unit/upper_log), unit)
        if p in examples:
            t = examples[p]
            assert total > threshold(t)
            print(f'Primes through {p}: reciprocal-log sum > {float(total):.6f}, '
                  f'threshold {threshold(t)}, exponent {t}.')


def digit_sum(z, p):
    total = 0
    while z:
        z, d = divmod(z, p)
        total += d
    return total


def check_carries():
    rows = 0
    for n in range(3, 81):
        scales = {}
        for p in PRIMES:
            m, length = 1, 0
            while m <= n:
                m *= p
                length += 1
            scales[p] = m, length
        for x in range(1, n//3+1):
            for y in range(x, (n-x)//2+1):
                z = n-x-y
                coefficient = comb(n, x)*comb(n-x, y)
                rows += 1
                for p, (m, length) in scales.items():
                    total = sum(digit_sum(v, p) for v in (x, y, z, m-1-n))
                    value, v = coefficient, 0
                    while value % p == 0:
                        value //= p
                        v += 1
                    assert total == (p-1)*(length+v)
    print(f'Four-word identities agree with binomial products on {rows} rows in five bases.')


if __name__ == '__main__':
    certify_assignments()
    certify_primorials()
    check_carries()
