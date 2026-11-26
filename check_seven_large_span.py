"""Check the constants and affine-line polynomial controls in the large-span theorem."""
from fractions import Fraction
from itertools import product
from math import gcd

from check_seven_refined_cutoffs import below_one
from check_translated_triple_lines import specialize


def constants():
    assert below_one([(2,1,101),(3,1,1),(2,300,1000),(11,102,1000),(11,-1,1)])
    assert 29 < 2**5 and 9**16 < 11**15
    # At H >= 10^14 these are rational upper bounds on normalized error terms.
    assert 40000**3 <= 10**14
    assert Fraction(13,40000) < Fraction(1,1800)
    assert Fraction(39,40000) + Fraction(22*3600*60,10**7) < Fraction(9,16)
    assert Fraction(39,40000) + Fraction(22*12*4,10**7) < Fraction(1,16)


def affine_controls():
    count = 0
    directions = [(i,j) for i in range(5) for j in range(-4,5) if gcd(i,j) == 1]
    for p, n, slopes, t, v in product((11,17), (3,9,15), ((2,4),(4,8),(8,10)),
                                      (1,3,7), (-8,-3,0,3,8)):
        h, k = slopes
        first = ((0,0),(n,0))
        second = ((0,0),(h,0),(k,0),(0,1))
        for direction in directions:
            f,g = (specialize(p, block, direction, (t,v)) for block in (first,second))
            assert f.gcd(g).degree() == 0
            count += 1
    print(f"{count} exact coprimality checks on exponent lines")


if __name__ == '__main__':
    constants()
    affine_controls()
