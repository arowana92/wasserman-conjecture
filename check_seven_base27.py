"""Rational bounds at 27 and product-of-roots controls for the index-three matching."""
from fractions import Fraction
from math import gcd

from check_gap_resultants import resultant
from check_seven_refined_cutoffs import below_one
from check_seven_large_bases import X


def bounds():
    b = 27
    rows = [
        [(3,3,2),(6,1,1),(9,1,1)],
        [(3,3,2),(6,1,1),(b+1,2,101),(9,1,1)],
        [(3,3,4),(6,1,2),(b+1,2,3),(9,1,1)],
        [(Fraction(11,3),9,8),(4,3,4),(b+1,1,2),(9,1,1)],
        [(32,1,3),(b+1,2,3),(9,1,1)],
        [(Fraction(19,4),2,5),(4,4,15),(b+1,14,15),(9,1,1)],
        [(27,1,1),(4,2,3),(b+1,2,3),(4,1,101),(b+1,4,101)],
        [(Fraction(27,8),1,101),(8,1,1),(10,1,1),(9,1,1)],
        [(Fraction(27,8),1,101),(8,1,1),(6,1,1),(b+1,2,101),(9,1,1)],
        [(5,3,4),(6,1,2),(b+1,2,3),(9,1,1)],
        [(48,1,1),(2,24,7),(Fraction(27,8),1,101),(2,22,707),(b+1,2,101)],
        [(3,3,2),(54,1,1),(b,18,101)],
        [(3,7,24),(54,1,1),(b,2,3)],
    ]
    for i, row in enumerate(rows):
        assert below_one(row+[(b,-2,1)]), i
    for n in (7,9,11,13,17,19):
        d, k = 5*n, n+4
        assert 9**101*32**(d-k)*(b+1)**(2*k) < b**202
    print(f"{len(rows)} uniform bounds and six short matching-index comparisons")


def root_products():
    count = 0
    for w in range(1,18):
        f = X**(3*w)+1
        for length in (1,3,5):
            for c in range(1,4*w+1):
                if gcd(c,w) != 1:
                    continue
                y = 1+X**(length*w)+X**c
                value = resultant(f,y)
                if length == 3:
                    assert value == 1
                else:
                    assert value <= 4*3**w
                count += 1
    print(f"{count} exact index-three resultant checks")


if __name__ == '__main__':
    bounds()
    root_products()
