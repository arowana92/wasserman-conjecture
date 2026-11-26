"""Check the constants and affine-line polynomial controls in the large-span theorem."""
from fractions import Fraction
from itertools import product
from math import gcd

from check_seven_refined_cutoffs import below_one
from check_translated_triple_lines import specialize


def constants():
    for row in ([(3,1,1),(10,1,2)], [(4,1,1),(6,1,2),(12,1,100000)],
                [(4,1,100000),(11,18,100000),(3,1,1),(2,2,3),(11,1,3)],
                [(3,1,1),(2,1,10),(11,7,15)]):
        assert below_one(row+[(11,-197,200)])
    assert 9**16 < 11**15
    assert 460**3 <= 10**8 and 250**2*1600 == 10**8
    assert Fraction(1,100)-Fraction(4,460) > Fraction(1,800)
    assert Fraction(12,460)+Fraction(8,250) < Fraction(1,16)


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
