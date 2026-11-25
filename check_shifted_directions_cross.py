"""Independent structural enumeration and exact arithmetic controls for shifted rows."""

from itertools import product
from math import comb, gcd

from carry_blocks import CarryGraph, intersection
from check_prime_power_cofactors import shapes
from check_shifted_directions import binomial_divisor, possible_scales
from check_shifted_prefixes import floor_cover
from search_endpoint import atoms, smallest_factors


def independent_shapes(m):
    # Enumerate the shifted multipliers first. Two zero remainders in C
    # give the other slot a uniquely determined residue. Construct C by CRT.
    factors = []
    rest = m
    for p in range(2, m+1):
        w = 1
        while rest % p == 0:
            w *= p
            rest //= p
        if w > 1:
            factors.append((p,w))
    answers = set()
    for a in range(1, (m-1)//2+1):
        for b in range(a, m-a):
            c = m-a-b
            if gcd(a,b,c) != 1:
                continue
            common = gcd(a,b)
            if common < 3:
                continue
            required = [w for p,w in factors if
                        max(gcd(a,w),gcd(b,w)) != w or min(gcd(a,p),gcd(b,p)) != 1]
            for u in range(3, common+1):
                if common % u or gcd(m,u) != 1 or any(2*w >= u for w in required):
                    continue
                for slots in product(range(3), repeat=len(required)):
                    residues = [0,0,0]
                    period = 1
                    for w,slot in zip(required,slots):
                        for i,r in enumerate(residues):
                            goal = u if i == slot else 0
                            residues[i] = next(x for x in range(r, period*w, period)
                                               if (x-goal) % w == 0)
                        period *= w
                    f0,g0,h0 = (r or period for r in residues)
                    for f in range(f0,u,period):
                        for g in range(g0,u-f,period):
                            h = u-f-g
                            if h >= h0 and (h-h0) % period == 0:
                                answers.add((u,a,b,c,*sorted((f,g,h))))
    return answers


def compare_enumerations():
    spf = smallest_factors(1000)
    total = surviving = 0
    # Include each of the two smallest surviving denominators and both
    # cases where allowing composite even scales introduces new shapes.
    for m in (*range(3,101), 364,380,738,986):
        factors = atoms(m,spf)
        original = {s for s in shapes(m,factors) if gcd(*s[1:4]) == 1}
        expected = independent_shapes(m)
        if original != expected:
            raise ArithmeticError(("shape mismatch",m,original ^ expected))
        total += len(expected)
        retained = {(s,parity) for s in expected
                    for parity, _, _ in possible_scales(m,factors,s,spf)}
        surviving += len(retained)
        if len(retained) != {364:1,380:1,738:1,986:3}.get(m,0):
            raise ArithmeticError(("local classification",m,retained))
    print(f"Reverse enumeration agrees on {total} primitive shapes and includes "
          f"{surviving} locally surviving shape/parity pairs.")


def exact_controls():
    ntests = accepted = shared = 0
    for m in range(3,17):
        spf = smallest_factors(m)
        for a in range(1, m):
            for b in range(a,m-a):
                c = m-a-b
                if gcd(a,b,c) != 1:
                    continue
                for scale in range(1,19):
                    parts = a*scale, b*scale, c*scale+1
                    n = sum(parts)
                    coefficient = comb(n,parts[0])*comb(n-parts[0],parts[1])
                    assert coefficient % binomial_divisor(a,b) == 0
                    for p,w in atoms(m,spf):
                        ntests += 1
                        if scale % p == 0:
                            shared += 1
                            assert coefficient % p == 0
                        if coefficient % p:
                            accepted += 1
                            assert (a+b)*scale % w == 1
                            assert (a % w == 0) != (b % w == 0)
    # Finite covers with positive results, computed directly by binomials.
    successes = 0
    for slopes in ((1,2,4),(1,3,2),(2,4,3)):
        for bases in ((2,), (3,5)):
            graphs = [CarryGraph(p,slopes,(0,0,1)) for p in bases]
            for mod, r in ((1,0),(3,1),(4,3)):
                direct = []
                for w in range(2,240):
                    if w % mod != r:
                        continue
                    a,b,c = slopes
                    n = sum(slopes)*w+1
                    k = comb(n,a*w)*comb(n-a*w,b*w)
                    if all(k % p for p in bases):
                        direct.append(w)
                assert intersection(graphs,239,mod,r)[0] == direct
                assert floor_cover(slopes,bases,mod,r,239)[0] == direct
                successes += len(direct)
    assert accepted and shared and successes
    print(f"{ntests} exact prime tests, including {shared} shared-prime scales and "
          f"{accepted} avoidances. Both finite covers recover {successes} accepted inputs.")


def check_power_residues():
    # Each entry gives the denominator and W's full residue class, as printed
    # by the sieve. Test each prime-power component separately by CRT.
    progressions = [
        (364,747,988),(380,823,2470),(589,5060,8246),(589,937,8246),
        (629,1016,8806),(629,5419,8806),(646,1711,4522),(703,3880,9842),
        (703,8801,9842),(710,4643,4970),(738,7250,9594),(780,919,1820),
        (850,5651,5950),(874,35,9614),(880,2801,2992),(946,531,6622),
        (950,4463,10450),(962,6349,10582),(986,4215,6902),(986,3501,6902),
        (986,10784,12818),
    ]
    spf = smallest_factors(max(mod for _, _, mod in progressions))
    for k, expected in ((2,set()), (3,{880,986}), (105,set())):
        possible = set()
        for m,r,mod in progressions:
            if all(any(pow(i,k,w) == r % w for i in range(w)) for _,w in atoms(mod,spf)):
                possible.add(m)
        assert possible == expected, (k,possible)
    print("Prime-power CRT confirms all 21 residue tests at exponents 2, 3 and 105.")


if __name__ == "__main__":
    compare_enumerations()
    exact_controls()
    check_power_residues()
