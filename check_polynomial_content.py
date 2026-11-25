"""Independently enumerate the polynomial splits occurring in the content obstruction."""

from itertools import product
from math import gcd
from sympy import Poly, symbols

X = symbols("X")


def split_checks():
    tests = mixed = proportional = 0
    for coefficients, m in (((1,1),7), ((1,2),7), ((2,1),5),
                            ((1,0,1),4), ((1,1,1),7), ((2,0,1),3)):
        H = Poly(sum(v*X**i for i,v in enumerate(coefficients)), X)
        assert H.is_irreducible
        G = m*H
        def nonnegative(f):
            return all(v >= 0 for v in f.all_coeffs())
        choices = [Poly(sum(v*X**i for i,v in enumerate(row)), X)
                   for row in product(*(range(v*m+1) for v in coefficients)) if any(row)]
        # Divide exact polynomials instead of generating the proposed form.
        whole = [f for f in choices if f.rem(H).is_zero]
        unshifted = [f for f in choices if (X*f-1).rem(H).is_zero]
        shifted = [f for f in choices if (X*f+1).rem(H).is_zero]
        for left,middle in product(whole, repeat=2):
            right = G-left-middle
            if right.is_zero or not nonnegative(right):
                continue
            assert right.rem(H).is_zero
            ks = [int(f.exquo(H).LC()) for f in (left,middle,right)]
            assert min(ks) > 0 and sum(ks) == m
            for P in (2,5,7):
                entries = [P*int(f.eval(P)) for f in (left,middle,right)]
                assert sum(entries)//gcd(*entries) == m//gcd(*ks)
            proportional += 1
        for left,right in product(unshifted,shifted):
            middle = G-left-right
            if middle not in whole:
                continue
            tests += 1
            b = int((X*left-1).exquo(H).nth(1))
            t = int(middle.exquo(H).nth(0))
            h = m-b-t
            assert coefficients[0] == 1
            D = (H-1).exquo(Poly(X,X))
            assert left == b*H-D and right == h*H+D
            for P in (2,3,4,5,7,11,13,16,19,25,31):
                value = int(H.eval(P))
                row = P*int(left.eval(P)), P*int(middle.eval(P)), P*int(right.eval(P))+1
                u = gcd(*row)
                if u < 3:
                    continue
                mixed += 1
                assert t % u == 0 and gcd(u,m*P*value) == 1
                assert b >= 2 and h >= 1 and m >= 7
                assert max(coefficients[i+1]-b*coefficients[i]
                           for i in range(len(coefficients)-1)) <= 0
                if m == 7:
                    power = 7
                    while (m*P*value) % (power*7) == 0:
                        power *= 7
                    if sorted(x % power for x in row) == [0,0,1]:
                        assert (h*P+1) % m == 0 and P % 7 in (3,6)
    assert tests and mixed and proportional
    print(f"Checked {proportional} proportional splits, {tests} mixed splits, "
          f"and {mixed} mixed specializations with row gcd at least three.")


if __name__ == "__main__":
    split_checks()
