"""Certify consecutive cubic weights and two quartic successors at every prime power."""

from collections import Counter
from itertools import combinations
from math import gcd, lcm

from flint import fmpz_poly

from check_close_powers import INITIAL_LIMIT
from check_polynomial_powers import (
    assignment_bounds, coeffs, evaluate, full_power_bound, power_limit, row_shapes, test_assignment,
)
from search_even_endpoint import prime_powers


def reduce_predecessor_bound(G, row, bound):
    """Remove the common polynomial scale before using two constant parts."""
    if len(bound) != 1:
        return bound
    polynomials = [fmpz_poly(list(part)) for part in row]
    shared = polynomials[0].gcd(polynomials[1]).gcd(polynomials[2])
    reduced = [f // shared for f in polynomials]
    constants = [int(f[0]) for f in reduced if f.degree() == 0]
    value = abs(bound[0])
    for a,b in combinations(constants,2):
        value = gcd(value,lcm(a,b))
        paired = (a+b)//gcd(a,b)
        for r,w in prime_powers(value).items():
            if paired % r == 0 or all((t*evaluate(G,t)-1) % r for t in range(r)):
                value //= w
    return (value,)


def verify(shifts, offsets):
    quotient = fmpz_poly([1])
    for h in shifts:
        quotient *= fmpz_poly([h,1])
    G = coeffs(quotient)
    lower = max(G)+1
    if (lower-1)*evaluate(G,lower-1)+1 > INITIAL_LIMIT:
        raise ArithmeticError("small-input interval needs a separate check")
    counts = Counter()
    for row in row_shapes(G):
        for eps in offsets:
            for shifted in ((0,1,2) if eps else (None,)):
                if shifted == 1 and row[0] == row[1]:
                    continue
                counts["polynomial rows"] += 1
                first,second = assignment_bounds(G,row,shifted)
                high = power_limit(G,eps,first,second)
                if high is None:
                    first = full_power_bound((1,)+G if eps else G,first)
                    second = full_power_bound(G if eps else (-1,)+G,second)
                    if not eps:
                        second = reduce_predecessor_bound(G,row,second)
                    high = power_limit(G,eps,first,second)
                if high is None:
                    raise ArithmeticError(("unbounded row",shifts,row,shifted,first,second))
                if high < lower or high*evaluate(G,high)+eps <= INITIAL_LIMIT:
                    continue
                for P in range(lower,high+1):
                    if P*evaluate(G,P)+eps > INITIAL_LIMIT:
                        counts["integer evaluations"] += 1
                        test_assignment(G,row,shifted,P,counts)
    print(f"Shifts {shifts}, offsets {offsets}: {dict(counts)}",flush=True)
    return counts


if __name__ == "__main__":
    total = Counter()
    for shifts in combinations(range(1,6),2):
        total.update(verify(shifts,(0,1)))
    for shifts in ((1,2,3),(1,2,4)):
        total.update(verify(shifts,(1,)))
    print("Totals:", dict(total))
