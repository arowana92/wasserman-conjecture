"""Enumerate every primitive shifted row at the divisors of 1651, 2222 and 2563."""

from math import gcd

from check_prime_power_cofactors import shapes
from check_shifted_directions import possible_scales
from check_shifted_directions_cross import independent_shapes
from search_endpoint import atoms, smallest_factors


def check():
    spf = smallest_factors(2563)
    all_cases = {}
    for m in range(3,2564):
        if 1651%m and 2222%m and 2563%m:
            continue
        fs = atoms(m,spf)
        forward = {shape for shape in shapes(m,fs) if gcd(*shape[1:4]) == 1}
        assert forward == independent_shapes(m)
        kept = []
        for shape in sorted(forward):
            for parity,progression,bases in possible_scales(m,fs,shape,spf):
                kept.append((shape,parity,progression,tuple(bases)))
        if kept:
            all_cases[m] = kept
        print('Denominator, primitive shapes, surviving parities:', m,len(forward),len(kept))
    assert all_cases == {
        1111: [((7,308,707,96,1,2,4),0,(7858,15554),(7,11,101)),
               ((7,308,707,96,1,2,4),1,(81,15554),(7,11,101))],
        1651: [((4,104,1524,23,1,1,2),1,(2297,6604),(2,13,127))],
        2222: [((7,77,2121,24,1,2,4),1,(8101,15554),(7,11,101))]}
    for s in (3,5):
        assert all(sum(a*t%s for a in (1,2,4)) >= s for t in range(1,s))
        for slopes in ((77,2121,24),(616,1414,192)):
            parts = slopes
            while sum(t%s for t in parts) < s:
                assert max(parts)
                parts = tuple(t//s for t in parts)


if __name__ == '__main__':
    check()
