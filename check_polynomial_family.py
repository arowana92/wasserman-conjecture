"""Apply the polynomial gcd criterion to P(2**k), with ascending coefficients."""

from argparse import ArgumentParser
from functools import lru_cache

from sympy import Poly, QQ, symbols


def check(coefficients):
    x = symbols("x")
    terms = [Poly((1 << bit)*x**degree, x, domain=QQ)
             for degree, coefficient in enumerate(coefficients)
             for bit in range(coefficient.bit_length()) if coefficient & (1 << bit)]
    total = sum(terms, Poly(0, x, domain=QQ))
    predecessor = total - 1
    full = (1 << len(terms)) - 1

    @lru_cache(None)
    def part(mask):
        if mask == 0:
            return Poly(0, x, domain=QQ)
        lowest = mask & -mask
        return part(mask-lowest) + terms[lowest.bit_length()-1]

    @lru_cache(None)
    def gcds(mask):
        poly = part(mask)
        return poly.gcd(total), poly.gcd(predecessor)

    partitions = tested = 0
    maximum = None
    witness = None
    for first in range(1, full, 2):  # The first part contains the lowest indexed term.
        rest = full ^ first
        second = (rest-1) & rest
        while second:
            third = rest ^ second
            if second < third:
                partitions += 1
                masks = (first, second, third)
                if all(part(mask).degree() > 0 for mask in masks):
                    tested += 1
                    (a, d), (b, e), (c, f) = map(gcds, masks)
                    at_total = a.gcd(b).gcd(c)
                    at_predecessor = d.gcd(e).lcm(d.gcd(f)).lcm(e.gcd(f))
                    degree = 3*at_total.degree() + 2*at_predecessor.degree()
                    if maximum is None or degree > maximum:
                        maximum = degree
                        witness = [part(mask).as_expr() for mask in masks]
            second = (second-1) & rest
    print(f"{total.as_expr()}: {partitions} partitions, {tested} with three nonconstant parts")
    if maximum is None:
        print("Every partition has a bounded part. The family is proved for sufficiently large k.")
    elif maximum < 2*total.degree():
        print(f"Maximum weighted gcd degree {maximum} < {2*total.degree()}. "
              "The family is proved for sufficiently large k.")
    else:
        print(f"Inconclusive at weighted gcd degree {maximum} for the split {witness}.")


if __name__ == "__main__":
    cli = ArgumentParser(description=__doc__)
    cli.add_argument("coefficients", nargs="+", type=int)
    coeff = cli.parse_args().coefficients
    if len(coeff) < 2 or coeff[-1] <= 0 or min(coeff) < 0:
        cli.error("use nonnegative integer coefficients and a positive leading coefficient")
    check(coeff)
