"""Power residues and carry columns for every remaining primitive denominator through 1000."""
from itertools import combinations
from math import isqrt, lcm
from fractions import Fraction

from check_prime_power_cofactors import shapes
from check_shifted_directions import possible_scales
from search_endpoint import atoms, smallest_factors

DENOMINATORS = (380,589,629,646,703,710,738,780,850,874,946,950,962,986)
PRIMES = tuple(p for p in range(3,80,2) if all(p%i for i in range(2,isqrt(p)+1)))


def residue_exclusion(R, M, exponent, spf):
    return any(R%w not in {pow(z,exponent,w) for z in range(w)} for _,w in atoms(M,spf))


def column_exclusion(slopes, R, p, h):
    # If W is a p^h-th power, its root is R mod p. Raising any lift of
    # that root to p^h fixes W modulo p^(h+1).
    power = p**(h+1)
    low = pow(R%p,p**h,power)
    for x in range(low,power*p,power):
        entries = (slopes[0]*x,slopes[1]*x,slopes[2]*x+1)
        free_by_layers = all(sum(z%(p**j) for z in entries) < p**j
                             for j in range(1,h+3))
        digits = list(entries)
        free_by_digits = True
        for _ in range(h+2):
            if sum(z%p for z in digits) >= p:
                free_by_digits = False
                break
            digits = [z//p for z in digits]
        assert free_by_digits == free_by_layers
        if free_by_layers:
            return False
    return True


def minimal_covers(sets):
    covers = {frozenset()}
    for obstruction in sets:
        candidates = {cover | {p} for cover in covers for p in obstruction}
        covers = {c for c in candidates if not any(b<c for b in candidates)}
    return covers


def verify():
    spf = smallest_factors(20000)
    records = []
    cube_shapes = []
    for m in DENOMINATORS:
        factors = atoms(m,spf)
        for shape in shapes(m,factors):
            for parity,(R,M),first in possible_scales(m,factors,shape,spf):
                assert M < len(spf)
                slopes = shape[1:4]
                excluded = {p for p in PRIMES if residue_exclusion(R,M,p,spf)
                            or (p in first and M%p == 0 and R%p != 0 and
                                column_exclusion(slopes,R,p,1))}
                if 3 not in excluded:
                    cube_shapes.append((m,slopes,R,M))
                    assert m == 986 and {7,17,29} <= set(first)
                    assert column_exclusion(slopes,R,29,2)
                records.append((parity,excluded))
    assert len(records) == 19
    assert sorted(cube_shapes) == [(986,(119,203,664),4215,6902),
                                    (986,(203,714,69),3501,6902)]
    total = minimal_covers([s for _,s in records])
    odd = minimal_covers([s for parity,s in records if parity])
    assert {tuple(sorted(c)) for c in total} == {
        (3,7),(3,17),(5,7,17),(7,11,13,17),(7,11,17,19,41),
        (7,13,17,19,37),(7,17,19,37,41)}
    assert {tuple(sorted(c)) for c in odd} == {
        (3,7),(3,17),(5,7,17),(7,11,13,17),(7,11,17,19),(7,17,19,37)}
    assert sorted({pow(z,7,29) for z in range(1,29)}) == [1,12,17,28]
    for m,slopes,R,M in cube_shapes:
        assert R%17 == 16 and R%29 in (10,21)
        print("Cube shape",slopes,"initial residue",R,"mod",M,
              "29th power excluded:",column_exclusion(slopes,R,29,1))
    guaranteed = (2,21,51,595,17017,24871,83657,2523)
    density = sum((Fraction((-1)**(k+1),lcm(*subset))
                   for k in range(1,len(guaranteed)+1)
                   for subset in combinations(guaranteed,k)), Fraction())
    assert density == Fraction(11490087461,21558946695)
    print("Density of simultaneous guaranteed exponents", density)
    print("Minimal prime covers for all primitive scales:",
          sorted(tuple(sorted(c)) for c in total))
    print("All 19 remaining shape/parity cases checked.")


if __name__ == '__main__':
    verify()
