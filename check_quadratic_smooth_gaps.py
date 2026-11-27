"""Check the heights and numerical margins in the quadratic smooth-gap theorem."""

from fractions import Fraction as F
from math import factorial, isqrt, prod

from check_seven_quadratic_bases import limits
from check_seven_quadratic_suffixes import RESIDUAL, upper_graph


def root_upper(n):
    unit = 10**6
    return F(isqrt(n*unit**2)+1, unit)


def matveev_upper(degree, heights, log_degree):
    t = len(heights)
    return F(7,5)*30**(t+3)*t**4*root_upper(t)*degree**2*(1+log_degree)*prod(heights)


def check():
    for D, _, _, z, w in RESIDUAL:
        _, graph = upper_graph(D, z, w)
        _, rays = limits(graph)
        for f in rays:
            for slope in (F(1), F(10,D), F(2), F(5)):
                c = F(8*D,5)*f*slope**2
                assert F(1,1000) < c < 10**6
                assert max(c.numerator, c.denominator) < 10**11
                # The 5-adic valuation of each rational is either -1 or +1.
                assert (c.numerator % 5 == 0) != (c.denominator % 5 == 0)
                assert c.numerator % 25 and c.denominator % 25

    # Certified natural-log upper bounds, using a positive Taylor lower sum for exp.
    upper = dict(zip((2,3,5,7,11,13,17,19,23),
                     (693148,1098613,1609438,1945911,2397896,2564950,2833214,2944439,3135495)))
    logs = {p:F(n,10**6) for p,n in upper.items()}
    for p, bound in logs.items():
        assert sum(bound**k / factorial(k) for k in range(35)) > p

    reference = 68*10**28
    for y, delta in ((17,F(1,10000)),(19,F(1,100)),(23,F(4,5))):
        heights = [F(26), *(a for p,a in logs.items() if p <= y)]
        coefficient = matveev_upper(1, heights, 0)
        assert delta*reference > 10 + 75*coefficient
        print(f'{y}-smooth: coefficient below {int(coefficient)+1}, delta={delta}.')

    for count, bits in ((1,10**13),(2,500000),(3,2000),(4,100)):
        heights = [F(26),logs[2]] + [bits*logs[2]]*count
        coefficient = matveev_upper(1, heights, 0)
        assert reference/100 > 10 + 75*coefficient
        print(f'{count} generators: bases through 2^{bits}.')

    # The roots of X^2-aX+e, e = +/-1, and both Lucas sequences.
    bits = 5*10**6
    root_log = (bits+2)*logs[2]
    heights = [2*logs[2], 2*(26 + 2*root_log), root_log]
    coefficient = matveev_upper(2, heights, logs[2])
    assert reference/100 > 10 + 75*coefficient
    print(f'Lucas recurrence coefficients through 2^{bits}.')


if __name__ == '__main__':
    check()
