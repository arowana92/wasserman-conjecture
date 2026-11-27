"""Independently read the middle binary columns from high to low."""

from fractions import Fraction
from functools import lru_cache
from math import lcm


DATA = (
    (4587, (3336, 308, 943), (Fraction(256,417),)),
    (7383, (1284, 552, 5547), (Fraction(19597,35328), Fraction(486233,876544),
                                Fraction(836139,1507328), Fraction(46678369,84148224))),
    (9453, (1644, 184, 7625), (Fraction(256,411),)),
)


def verify_shape(D, slopes, rays):
    a, b, c = slopes
    assert sum(slopes) == D
    multiple = lcm(a, b, c)
    dx, dy, dz = (multiple//s for s in slopes)

    def possible(x, y, z):
        return max(0, x*dx, y*dy, z*dz) < min(multiple, (x+1)*dx, (y+1)*dy, (z+2)*dz)

    stops = set()
    for x in range(a):
        for y in range((x*b)//a, ((x+1)*b)//a+1):
            for z in range((x*c)//a-1, ((x+1)*c)//a+1):
                if possible(x,y,z) and z >= 0 and not ((x&y) | (x&z) | (y&z)):
                    stops.add((x,y,z))

    @lru_cache(maxsize=None)
    def lower(state, digit, shift):
        x = 2*state[0] - a*digit
        y = 2*state[1] - b*digit
        z = 2*state[2] - c*digit + shift
        choices = ((x,y,z), (x+1,y,z), (x,y+1,z), (x,y,z+1))
        return tuple(s for s in choices if possible(*s))

    results = {}
    for f in rays:
        for parity in (0, 1):
            radicand = Fraction(8*D*(2**parity),5) * f
            # Recover a root prefix by testing each possible next binary digit.
            depth = 256
            bound = 100 * radicand * 2**(2*depth)
            high = 1
            while 4**high <= bound:
                high += 1
            root = 0
            for i in reversed(range(high)):
                trial = root + (1 << i)
                if trial*trial < bound:
                    root = trial
            assert (root+Fraction(1,2**400))**2 < bound < (root+1-Fraction(1,2**400))**2
            order = 0
            while 4**(order+1) < radicand:
                order += 1
            assert 4**order < radicand < 4**(order+1)
            shift_order = high-depth-1

            length = 40
            prefix = int(f*2**length)
            assert length == prefix.bit_length()
            for endpoint in (prefix, prefix+1):
                assert abs(f*2**length-endpoint) * 2**(100-length) > 210
            states = frozenset(stops)
            for digit in bin(prefix)[2:]:
                states = frozenset(child for s in states for child in lower(s,int(digit),0))
            first = states
            positions = {}
            period = 0
            while True:
                for s in states:
                    positions[s] = positions.get(s,0) | (1 << period)
                period += 1
                states = frozenset(child for s in states for digit in (0,1)
                                   for child in lower(s,digit,0))
                if states == first:
                    break
                assert period < 1200
            for digit in bin(root)[2:]:
                after = {}
                for s, mask in positions.items():
                    for x in (0,1):
                        for child in lower(s,x,int(digit)):
                            after[child] = after.get(child,0) | mask
                positions = after
            mask = 0
            for val in positions.values():
                mask |= val
            allowed = {(k+length+shift_order+order-parity-4) % period
                       for k in range(period) if mask >> k & 1}
            results[D,f,parity] = period, allowed
    return results


if __name__ == '__main__':
    results = {}
    for D, slopes, rays in DATA:
        results.update(verify_shape(D, slopes, rays))
    from check_seven_quadratic_middle import check
    assert results == check()
    print('All twelve full sets of order classes agree.')
