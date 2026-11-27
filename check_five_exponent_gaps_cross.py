"""Verify exponent-gap certificates by Lucas digits and factorial floor intervals."""
import json
import sys
from fractions import Fraction
from heapq import merge
from math import comb, factorial
from pathlib import Path


def lucas(t):
    top, bottom = 5*t + 1, 3*t
    while bottom:
        if bottom % 5 > top % 5:
            return False
        top //= 5
        bottom //= 5
    return True


def lucas_prefixes(multiplier, depth):
    # The preceding input digit supplies this column of 5T+1. Carry the
    # products 3T and multiplier*T from below, consuming one digit of T.
    states = [(0, 1, 1, 0, 0, 0)]
    place = 1
    for column in range(depth):
        following = []
        for t, prev, prev_y, triple, triple_y, product in states:
            for digit in ((2,) if column == 0 else range(5)):
                high, y = divmod(product + digit*multiplier, 5)
                carry, a = divmod(3*digit + triple, 5)
                carry_y, b = divmod(3*y + triple_y, 5)
                if a <= prev and b <= prev_y:
                    following.append((t + place*digit, digit, y, carry, carry_y, high))
        states = following
        if not states:
            break
        place *= 5
    return [state[0] for state in states]


def floor_cells(front, unit, slack):
    # In a cell floor(5y)=k, find the least y that can possibly satisfy
    # floor(3y)+floor(2y+epsilon)>=k, then include the closed upper end.
    following = []
    width = 6*unit
    for left, right in front:
        for k in range(left//width - 1, right//width + 1):
            j = (15*unit*k - slack) // (25*unit)
            threshold = min(max(10*unit*a, 15*unit*(k-a) - slack)
                            for a in (j, j+1))
            lo = max(left, width*k, threshold)
            hi = min(right, width*(k+1))
            if lo <= hi:
                if following and lo <= following[-1][1]:
                    following[-1] = (following[-1][0], max(hi, following[-1][1]))
                else:
                    following.append((lo, hi))
    return following


def large_exclusion(multiplier, height):
    precision = 5**(height+7)
    scale = 30*multiplier*precision
    near = 1
    while 5*near <= multiplier:
        near *= 5
    intervals = [(scale//5, scale)]
    place = 5*precision
    for _ in range(height+2):
        intervals = floor_cells(intervals, multiplier*place, 15*multiplier)
        intervals = floor_cells(intervals, near*place, 15)
        place //= 5
        if not intervals:
            return True
    return False


def density_checks():
    # Multiply surviving unit residues directly, without discrete logarithms.
    mod = 5**9
    units = lucas_prefixes(1, 9)
    ratios = {t*inverse % mod for u in units for inverse in (pow(u,-1,mod),) for t in units}
    assert len(units) == 1145 and len(ratios) == 151047
    # Construct the leading union from floor cells. Stream quotient intervals
    # in increasing order using a heap of sorted lists; all endpoints are exact.
    scale = 30*5**6
    cells = [(scale//5, scale)]
    for j in range(7):
        cells = floor_cells(cells, 5**(6-j), 0)
    def quotients(c, d, m):
        for a, b in cells:
            if m*b >= c and m*a <= 5*d:
                yield max(Fraction(1), Fraction(m*a,d)), min(Fraction(5), Fraction(m*b,c))
    grid = 10**12
    total = 0
    right = Fraction(1)
    for lo, hi in merge(*(quotients(c,d,m) for c,d in cells for m in (1,5))):
        if lo > right:
            increment = grid*(lo-right)/lo
            total += increment.numerator // increment.denominator
        right = max(right, hi)
    assert right == 5 and total*10000 > 39*161*grid
    assert sum(Fraction(161,100)**j/factorial(j) for j in range(10)) > 5
    assert Fraction(151047,4*5**8)*Fraction(61,100) < Fraction(6,100)
    print('Independent quotient and ratio bounds verified.')


def main(path):
    proofs = {int(k):v for k,v in json.loads(path.read_text()).items()}
    assert set(proofs) == set(range(4,4001,4))
    for t in range(1,160):
        assert lucas(t) == (comb(5*t+1, 3*t) % 5 != 0)
    # An actual simultaneous success for a different multiplier checks both
    # digit and interval logic on a nonempty example.
    control = next((t,d) for t in range(2,100) for d in range(6,100,5)
                   if lucas(t) and lucas(t*d))
    t, d = control
    assert t in lucas_prefixes(d, 5)
    large_t = (5**11 + 1) // 3
    assert large_t > 5**(3+7) and lucas(large_t) and lucas(26*large_t)
    assert not large_exclusion(26, 3)
    open_gaps = []
    for gap, record in proofs.items():
        multiplier = 1 << gap
        if record[0] == 'open':
            open_gaps.append(gap)
        elif record[0] == 'low':
            assert not lucas_prefixes(multiplier, record[1])
        else:
            assert record[0] == 'high'
            height, count = record[1:]
            assert large_exclusion(multiplier, height), (gap,height)
            short = lucas_prefixes(multiplier, height+7)
            assert len(short) == count
            assert not any(lucas(t) and lucas(multiplier*t) for t in short)
    assert set(open_gaps) == {1356,1424,2936,3364}
    print('All 1,000 records checked. Positive Lucas control:', control)
    density_checks()


if __name__ == '__main__':
    main(Path(sys.argv[1]))
