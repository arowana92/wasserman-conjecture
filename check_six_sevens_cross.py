"""Independent six-sevens certificate with standard-library integer intervals."""
from bisect import bisect_left
from math import gcd, isqrt


class ScaledLogs:
    """Enclose scale*log(x) using only integer arithmetic and an atanh series."""
    def __init__(self, bits):
        self.bits = bits
        self.scale = 1 << bits
        self.table = [self.series(i, 256) for i in range(256, 513)]
        self.two = self.table[-1]
        self.seven = self.log(7)

    def series(self, n, d):
        # 1 <= n/d <= 2, hence t <= 1/3. The integer terms below
        # underestimate scale*t^(2j+1) by at most 3 at every step.
        assert d <= n <= 2*d
        a, b = n-d, n+d
        if a == 0:
            return 0, 0
        term = (a << self.bits)//b
        square = (a*a << self.bits)//(b*b)
        total = steps = 0
        while term:
            total += 2*(term//(2*steps+1))
            term = term*square >> self.bits
            steps += 1
        return total, total+8*(steps+1)

    def log(self, n):
        k = n.bit_length()-1
        base = (n*256) >> k
        lo, hi = self.series(n*256, base << k)
        left, right = self.table[base-256]
        return k*self.two[0]+left+lo, k*self.two[1]+right+hi

    def base_seven(self, n):
        a, b = self.log(n)
        lo, hi = self.seven
        return (a << self.bits)//hi, -(-(b << self.bits)//lo)


def cover(bound, depth, threshold):
    logs = ScaledLogs(2*bound.bit_length()+192)
    scale = logs.scale
    left, right = logs.base_seven(2)
    # Continued fractions of the lower endpoint, followed by exact inequalities.
    n, d = left, scale
    h0, h1, q0, q1 = 0, 1, 1, 0
    choices = []
    while d:
        t, rem = divmod(n, d)
        n, d = d, rem
        h0, h1 = h1, h0+t*h1
        q0, q1 = q1, q0+t*q1
        if q1 > bound*10**27:
            break
        err = max(abs(q1*left-h1*scale), abs(q1*right-h1*scale))
        if q1 > 1000*(bound+3000) and err*64*(bound+3000) < scale:
            assert gcd(q1, h1) == 1
            choices.append(q1)
    assert choices
    powers = [1]
    while powers[-1] < choices[-1]*12800:
        powers.append(7*powers[-1])

    current = [(1, 0)]
    checked = largest_gap = 0
    for j in range(1, depth+1):
        following = []
        largest_gap = 0
        for coefficient, last in current:
            checked += 1
            odd = coefficient
            while odd % 7 == 0:
                odd //= 7
            if odd & (odd-1) == 0:
                assert odd.bit_length()-1 <= threshold
                q = choices[0]
            else:
                x, y = logs.base_seven(coefficient)
                for q in choices:
                    lower, upper = q*x, q*y
                    preceding = (lower//scale)*scale
                    if min(lower-preceding, preceding+scale-upper)*32 > scale:
                        break
                else:
                    raise ValueError(("unseparated logarithm", coefficient))
            assert q*12800 < 1 << (threshold//2)
            limit = bisect_left(powers, q*12800)
            largest_gap = max(largest_gap, limit)
            if j < depth:
                for gap in range(last, limit):
                    following.append((coefficient*powers[gap-last]+1, gap))
        current = following
    return largest_gap, checked


def digit_sum(n):
    s = 0
    while n:
        n, d = divmod(n, 7)
        s += d
    return s


def first_sum_at_least(lower):
    # Divide by seven first since every allowed summand is a positive power.
    n = max(1, (lower+6)//7)
    place = 1
    while True:
        rounded = ((n+place-1)//place)*place
        if digit_sum(rounded) <= 6:
            return 7*rounded
        place *= 7


def main():
    bounds = [10**180, 10**74, 10**55, 10**38, 10**21, 3000, 200]
    total = 0
    for depth, bound, target in zip((2, 3, 4, 5, 6, 6), bounds, bounds[1:]):
        threshold = min(target, 3000)
        last, count = cover(bound, depth, threshold)
        if depth < 6:
            c = 10**14*(2+3*(len(str(bound))-1))
            assert (last+2)*c**(6-depth) < target
        total += count
        print(f"Integer enclosure pass {depth} at {bound}: {count:,} checks.", flush=True)
    for e in range(17, 201):
        n = 1 << e
        radius = isqrt(16*n)
        least = first_sum_at_least(n-radius)
        if e == 17:
            assert least == 5*7**4+7**6 == 129654 and least <= n+radius
        else:
            assert least > n+radius
    possible = []
    for m in range(4):
        for k in range(6):
            if 2*k > 3*m+2 and 6+3*m+2*k <= 17:
                t = (32*8**m*(14*4**k+1)-18)//42
                assert t % 1008 != 587
                possible.append(t)
    assert possible == [171, 683, 2731, 10923, 5467, 21851]
    print(f"All {total:,} checks passed. All six terminal pair candidates excluded.")


if __name__ == "__main__":
    main()
