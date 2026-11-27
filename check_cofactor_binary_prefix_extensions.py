"""Check finite binary windows and a simultaneous class at 32 cofactors."""

from math import gcd, isqrt, lcm


def check_words():
    count = 0
    for a in range(1, 11):
        for b in range(a, 12):
            if not (a | b) & 1:
                continue
            h = max((w & -w).bit_length()-1 for w in (a, b, a+b))
            for first in (0, 1, 6, 17):
                for second in (0, 3, 10):
                    for length in range(1, 8):
                        mod = 1 << length
                        words = sum(not ((a*x+first) & (b*x+second) & (mod-1)) for x in range(mod))
                        # floor(phi^n) from phi^n = F_n phi + F_(n-1).
                        f, g = 0, 1
                        for _ in range(length+h):
                            f, g = g, f+g
                        bound = g-f + (f+isqrt(5*f*f))//2
                        assert words <= bound
                        count += 1
    return count


def overlap(x, y, bits, required):
    common = x & y & ((1 << bits)-1)
    for _ in range(required):
        if not common:
            return False
        low = common & -common
        layer = 2*low
        assert (x+y)//layer > x//layer + y//layer
        common -= low
    return True


def explicit_class(r, bits, max_d, weight_bound, required):
    # Test every grouped cut, including the unit offset on its complementary part.
    for a in range(1, max_d):
        for b in range(1, max_d-a+1):
            assert overlap(a*r, b*r+1, bits, required)
    count = even = 0
    for d in range(1, max_d+1):
        n = d*r+1
        for m in range(3, weight_bound*d+1):
            low = m & -m
            if gcd(m, d) != 1 or n % low:
                continue
            precision = bits - (low.bit_length()-1)
            mod = 1 << precision
            v = (n//low * pow(m//low, -1, mod)) % mod
            assert (m*v-n) % (1 << bits) == 0
            # All pairs of normalized positive parts that leave a positive third.
            for a in range(1, (m-1)//2+1):
                for b in range(a, m-a):
                    assert overlap(a*v, b*v, precision, required)
                    count += 1
                    even += (m%2 == 0)
    return count, even


def check_lifts():
    count = 0
    for d in range(1, 13):
        for m in range(3, 21):
            if gcd(d, m) != 1:
                continue
            a = (m & -m).bit_length()-1
            g = (d & -d).bit_length()-1
            t, length = 5, 5
            for r in range(1, 1 << t, 2):
                if (d*r+1) % (1 << a):
                    continue
                start = t+g-a
                mod = 1 << (start+length)
                uppers, lowers = set(), set()
                for z in range(1 << length):
                    q = r + (z << t)
                    v = ((d*q+1) >> a) * pow(m >> a, -1, mod) % mod
                    uppers.add(v >> start)
                    lowers.add(v % (1 << start))
                assert len(lowers) == 1 and len(uppers) == 1 << length
                count += 1
    return count


def direct_lifts():
    counts = [0, 0]
    for p in (2, 3, 5):
        for d in range(2, 15):
            bound = 2*d
            L = lcm(*(s for s in range(1, bound+1) if gcd(s, p*d) == 1))
            g = 0
            while p**(g+1) <= bound:
                g += 1
            t = 2
            for r in range(1, p**t):
                if r % p == 0:
                    continue
                v = -r % p**t
                v += p**t * ((pow(d, -1, L)-v)*pow(p**t, -1, L) % L)
                if v < 2:
                    v += p**t * L
                j = 1
                while p**j <= d*v:
                    j += 1
                H = max(t, j+1+g)
                for Q in (-v % p**H, p**H+(-v % p**H)):
                    assert Q % p**t == r
                    N = d*Q+1
                    for layer in (p**j, p**(j+1)):
                        for a in range(1, d):
                            assert N//layer - (a*Q)//layer - ((d-a)*Q+1)//layer == 1
                            counts[0] += 1
                        for m in range(3, bound+1):
                            if N % m:
                                continue
                            V = N//m
                            for a in range(1, m//3+1):
                                for b in range(a, (m-a)//2+1):
                                    parts = (a*V, b*V, (m-a-b)*V)
                                    assert N//layer - sum(x//layer for x in parts) == 2
                                    counts[1] += 1
    return counts


if __name__ == "__main__":
    print("Finite word counts:", check_words())
    print("Bijections after division by the row denominator:", check_lifts())
    counts = explicit_class(1058742659284095799, 64, 32, 2, 3)
    assert counts[0] == 94161 and counts[1] > 0
    print("Companion pairs checked, including those at even denominators:", counts)
    print("Direct prime-power layers in cuts and companion rows:", direct_lifts())
