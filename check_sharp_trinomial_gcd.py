"""Exact controls for the two-thirds exponent and its exceptional binary pair."""
from fractions import Fraction
from itertools import product
from math import gcd
from random import Random

from flint import fmpq_poly
from z3 import And, Not, Or, Real, Solver, unsat


EPS = Fraction(1, 10000)
X = fmpq_poly([0, 1])


def geometry():
    # Independent real-linear verification. A gap is shorter than T iff
    # every consecutive pair is closer than T. Every nonmaximal occupied
    # point therefore needs a successor within T.
    a, b, d = (Real(k) for k in ("a", "b", "d"))
    one = Fraction(1)
    third = one/3
    cut = third+EPS

    def short_gaps(points):
        return And(*(Or(x == 1, *(And(y > x, y-x <= cut) for y in points))
                     for x in points))

    solver = Solver()
    solver.add(a >= 0, a <= 1, b >= 0, b <= d, d <= 1, d >= 2*third-EPS,
               short_gaps([one*0, a, one, b, d]),
               short_gaps([one*0, 1-a, one, d-b, d]))

    def near(x,y,error):
        return And(x-y <= error, y-x <= error)

    def edge(mid, small_mid):
        return And(mid >= third-2*EPS, mid <= cut,
                   near(1-3*mid, 0, 6*EPS), near(d-small_mid-mid, 0, 3*EPS),
                   small_mid >= mid, small_mid <= 2*cut)

    middle = And(b >= third-5*EPS, b <= third+4*EPS,
                 near(d-2*b, 0, 9*EPS), near(1-3*b, 0, 15*EPS))
    solver.push()
    solver.add(d > 2*cut, Not(Or(edge(a,b), edge(1-a,d-b))))
    assert solver.check() == unsat
    solver.pop()
    solver.add(d <= 2*cut, Not(middle))
    assert solver.check() == unsat
    print("Both geometric implications over the reals verified with rational linear arithmetic")


def thresholds():
    third = Fraction(1,3)
    beta = 2*third-EPS
    q, height = 16, 10000
    edge_w, middle_w = third+EPS, third+4*EPS
    assert (Fraction(3,2)*edge_w+33*q*EPS+3*edge_w/q+
            Fraction(16*q+2,height) < beta)
    assert middle_w+66*q*EPS+2*middle_w/q+Fraction(16*q+2,height) < beta
    assert (2*(third+EPS)/(third-2*EPS)+Fraction(1,q) < 3 and
            1/(third-5*EPS)+Fraction(1,q) < 4)
    h = 10**24
    assert h == (10**6)**4 == (10**8)**3
    assert Fraction(1,100000)-Fraction(6,10**8)-Fraction(2,h) > Fraction(1,101000)
    losses = Fraction(5,200000)+Fraction(8,10**6)+Fraction(5,h)
    assert EPS-losses > Fraction(1,20000)
    assert third-losses > Fraction(11,35)
    assert Fraction(2*32+2,3*32+2) == Fraction(33,49)
    print("Height, exponent and matching thresholds verified")


def line_controls():
    bases = (2,3,5,7)
    count = exceptional = 0
    for p in bases:
        coeff = [Fraction(p)**j for j in range(-3,4)]
        for b, d, e in product(coeff, repeat=3):
            # Convert rationals through numerator and denominator for FLINT.
            f = 1+X+(X**3*b.numerator)/b.denominator
            quad = 1+X+(X**2*b.numerator)/b.denominator
            for u in range(4):
                extra = (X**u*d.numerator)/d.denominator
                cubic = 1+extra+(X**3*e.numerator)/e.denominator
                assert quad.gcd(cubic).degree() <= 1
                count += 1
                if u not in (1,2):
                    continue
                companion = 1+extra+(extra*X*e.numerator)/e.denominator
                common = f.gcd(companion)
                if p == 2 and (b,d,e,u) == (4,1,2,2):
                    assert common.degree() == 2
                    exceptional += 1
                else:
                    assert common.degree() <= 1
                count += 1
    assert exceptional == 1
    print(f"{count} exact small-slope polynomial gcds, one quadratic exception")


def numerical_controls():
    count = 0
    for w in (1,2,3,7,43,501,5001,20000):
        z = 1 << w
        assert gcd(1+z+4*z**3, 1+z**2+2*z**3) == 1-z+2*z*z
    rng = Random(785)
    for p in (2,3,5,7):
        for step in (40,500,10000):
            for _ in range(75):
                mid = rng.randint(step,2*step)
                off = rng.randint(-2,2)
                big = (0,step,3*step+off)
                for small in ((0,mid,mid+step+rng.randint(-2,2)),
                              (0,step,2*step+rng.randint(-2,2))):
                    h = max(big[-1],small[-1])
                    if {big, small} == {(0,step,3*step+2),(0,2*step,3*step+1)} and p==2:
                        continue
                    for first, second in ((big,small),
                                          (tuple(big[-1]-e for e in big),
                                           tuple(small[-1]-e for e in small))):
                        value = gcd(sum(p**e for e in first),sum(p**e for e in second))
                        # Sufficient integral rounding of the claimed bound.
                        if value > 4*p**((19997*h)//30000+int(h < 10000)):
                            raise AssertionError((p,first,second,value))
                        count += 1
    print(f"Eight exact members of the binary family and {count} nearby integer gcds")


if __name__ == '__main__':
    geometry()
    thresholds()
    line_controls()
    numerical_controls()
