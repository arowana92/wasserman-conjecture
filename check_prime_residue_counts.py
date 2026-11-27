"""Rational entropy certificates for mixed-base exceptional-weight bounds."""

from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import prod

from check_binary_windows import ln_bounds


# Prime, followed by the eight rational parameters in units of 1/10000.
PARAMETERS = (
    (3, 1158, 6570, 7793, 4443, 2305, 4309, 4030, 5020),
    (5, 1712, 6076, 6955, 5887, 2447, 5607, 3830, 6118),
    (7, 2123, 5689, 6425, 6802, 2546, 6410, 3680, 6785),
    (11, 2704, 5111, 5802, 7864, 2673, 7348, 3470, 7572),
    (13, 2923, 4886, 5603, 8193, 2717, 7647, 3390, 7828),
    (17, 3274, 4514, 5318, 8644, 2781, 8072, 3270, 8193),
    (19, 3419, 4358, 5213, 8804, 2806, 8228, 3220, 8331),
    (23, 3666, 4087, 5048, 9043, 2846, 8472, 3130, 8550),
    (29, 3962, 3756, 4873, 9280, 2889, 8727, 3030, 8781),
    (31, 4046, 3661, 4828, 9338, 2900, 8793, 3000, 8843),
    (37, 4266, 3409, 4718, 9474, 2928, 8954, 2930, 8989),
    (41, 4392, 3263, 4660, 9541, 2943, 9038, 2890, 9067),
    (43, 4449, 3196, 4635, 9570, 2949, 9074, 2870, 9102),
    (47, 4556, 3070, 4590, 9619, 2961, 9139, 2840, 9161),
    (53, 4699, 2901, 4534, 9678, 2976, 9220, 2790, 9239),
    (59, 4824, 2751, 4489, 9724, 2988, 9286, 2750, 9302),
    (61, 4863, 2705, 4476, 9737, 2992, 9306, 2740, 9320),
    (67, 4970, 2574, 4440, 9771, 3001, 9358, 2700, 9372),
    (71, 5036, 2494, 4419, 9790, 3007, 9388, 2680, 9401),
    (73, 5067, 2456, 4410, 9798, 3009, 9402, 2670, 9415),
    (79, 5155, 2347, 4384, 9821, 3016, 9440, 2640, 9452),
)

# Outward dyadic rounding before summing the logarithm series keeps the operands small.
@lru_cache(maxsize=None)
def log_interval(x):
    assert x >= 1
    two_power = max(0, x.numerator.bit_length() - x.denominator.bit_length())
    if 2**two_power > x:
        two_power -= 1
    m = x/(2**two_power)
    assert 1 <= m < 2
    unit = 1 << 50
    n = int(m*unit)
    low, _ = ln_bounds(F(n,unit),20)
    _, high = ln_bounds(F(n+1,unit),20)
    log2_lo, log2_hi = ln_bounds(F(2),20)
    return low + two_power*log2_lo, high + two_power*log2_hi


def log_up(x, base):
    return log_interval(x)[1]/log_interval(F(base))[0]


def costs(p, u, v):
    d = log_up(1+u,2) + log_up((1-v**p)/(1-v),p)
    return max(F(0),d-1), log_up(1/u,2), (p-1)*log_up(1/v,p)


def bounds(row):
    p, *numbers = row
    a,b,u,v,x,y,U,V = (F(n,10000) for n in numbers)
    c = 1-a-b
    assert min(a,b,c) > 0
    delta, alpha, beta = costs(p,u,v)
    base = a*log_up(1+x,2) + b*log_up((1-y**p)/(1-y),p)
    worst = max(base + c*(delta + (1-r)*alpha + (1-s)*beta)
                + (a*r/3)*log_up(1/x,2) + (b*(p-1)*s/3)*log_up(1/y,p)
                for r,s in product((0,1), repeat=2))
    separate = (worst+2*a/3)/(1+2*a/3)
    delta, alpha, beta = costs(p,U,V)
    together = (1+2*(delta+(alpha+beta)/4))/3
    return separate,together


def threshold(p):
    for cutoff, upper in ((3,8178), (5,8293), (7,8369), (11,8462),
                          (13,8493), (23,8581), (53,8665), (73,8685)):
        if p <= cutoff:
            return F(upper,10000)
    return F(86885,100000)


def check():
    primes = [p for p in range(3,80) if all(p%k for k in range(2,p))]
    assert primes == [row[0] for row in PARAMETERS]
    for row in PARAMETERS:
        p = row[0]
        separate,together = bounds(row)
        target = threshold(p)
        assert max(F(8113,10000), separate, together) < target, p
        if p == 3:
            assert together < F(7611,10000)
        print(p, f'{float(separate):.7f}', f'{float(together):.7f}', target)
    density = 1-prod(F(p-2,p) for p in primes)
    assert density > F(959,1000)
    print(f'Covered natural density {float(density):.9f}.')


if __name__ == '__main__':
    check()
