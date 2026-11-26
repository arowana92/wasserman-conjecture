"""Uniform finite counts of binary orders by simultaneous rotation separation."""
from functools import lru_cache
from itertools import product
from math import isqrt

from flint import fmpz_mat

from carry_blocks import CarryGraph
from check_logarithmic_intervals import (
    ceiling, exponent_floor, log_integer, log_ratio_bounds, prefix_counts,
)


@lru_cache(maxsize=None)
def separation(p, r, H):
    """Certify that the first H orbit points are separated by 1/D in sup norm."""
    # T is a large common denominator for the two rotation coordinates.
    T = H**4
    centers = []
    for base in (p,r):
        lo, hi = log_ratio_bounds(2, base, T.bit_length()+100)
        mid = T*(lo+hi)/2
        center = mid.numerator//mid.denominator
        assert center-1 <= T*lo <= T*hi <= center+2
        centers.append(center)

    D = isqrt(H)
    while True:
        # |n| < H and ||n*theta_i|| <= 1/D would give a short lattice vector.
        weight = ceiling(T,H*D)
        original = fmpz_mat([[weight,*centers],[0,T,0],[0,0,T]])
        reduced, change = original.lll(transform=True)
        assert abs(change.det()) == 1 and change*original == reduced
        u,v,_ = [[int(x) for x in row] for row in reduced.tolist()]
        cross = [u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0]]
        norm = sum(x*x for x in u)
        area = sum(x*x for x in cross)
        # Squared Gram-Schmidt lengths of the three basis vectors.
        lengths = ((norm,1), (area,norm), (int(reduced.det())**2,area))
        threshold = (H*weight)**2 + 2*(ceiling(T,D)+2*H)**2
        if all(num > threshold*den for num,den in lengths):
            return D
        D *= 2


@lru_cache(maxsize=None)
def counts(base, slopes, length):
    return prefix_counts(CarryGraph(base,slopes,(0,1)),length)


def box_bound(H, p, r, counts_for):
    D = separation(p,r,H)
    bases = (2,p,r)
    ranges = [range(max(1,exponent_floor(D,b)-3),exponent_floor(D,b)+5) for b in bases]
    prefixes = [fn(max(rg)) for fn,rg in zip(counts_for,ranges)]
    log_lower = {b: log_integer(b,40)[0] for b in (p,r)}
    best = H
    for lengths in product(*ranges):
        P,Q,R = (b**(n-1) for b,n in zip(bases,lengths))
        blocks = 1
        for b,other in ((p,Q),(r,R)):
            # floor(length * D) + 1 subarcs of length strictly less than 1/D.
            blocks *= 1 + ((D*(P+other) << 40) // (P*other*log_lower[b]))
        words = 1
        for table,n in zip(prefixes,lengths):
            words *= table[n]
        small = max(P,Q,R).bit_length()
        best = min(best,small+words*blocks)
    return best,D


def exact_separation_controls():
    """Check the rotation distances using just ordinary integer powers."""
    for H in (25, 100):
        for p, r in ((3, 5), (3, 7), (5, 7)):
            D = separation(p, r, H)
            powers = {b: 1 for b in (p, r)}
            for n in range(1, H):
                two = 1 << n
                target = two**D
                close = []
                for b in (p, r):
                    while powers[b]*b <= two:
                        powers[b] *= b
                    close.append(any(x**D <= b*target and target <= b*x**D
                                     for x in (powers[b], b*powers[b])))
                assert not all(close)
    print("All six finite separation controls passed by direct integer powers.")


def verify():
    from math import comb
    from check_logarithmic_intervals import forced_pairs, row_bound
    exact_separation_controls()
    cuts = ((5,3,3,5), (7,3,3,7), (7,5,5,7), (7,6,3,7),
            (9,5,3,5), (9,7,3,7))
    for exponent in (30,100,300,1000):
        H = 10**exponent
        totals = {}
        for d,k,p,r in cuts:
            slopes = (k,d-k)
            funcs = [lambda n,b=b: counts(b,slopes,n) for b in (2,p,r)]
            if (d,k) == (7,6):
                funcs[2] = lambda n: [comb(i+5,5) for i in range(n+1)]
            bound, D = box_bound(H,p,r,funcs)
            totals[d] = totals.get(d,0)+bound
            print(f"H=10^{exponent} case ({d},{k}), separation 1/({D//isqrt(H)} sqrt(H)),"
                  f" count <= {bound}",flush=True)
        total = min(H,sum(totals.values()))
        print(f"Combined bound for 5, 7 and odd 9 weights: {total}",flush=True)
        if exponent in (100, 300, 1000):
            ceiling_bound = {100: 232*10**97, 300: 160*10**291, 1000: 502*10**969}[exponent]
            assert total < ceiling_bound
        if exponent == 1000:
            bounded = sum(((p-1).bit_length()+1)*row_bound(sl,p,r,H+2)
                          for _,sl,p,r in forced_pairs())
            assert bounded < 145*10**778
            assert (H-total-bounded) // (total+bounded+1) >= 10**28
            print("Together with every prime-power cofactor through 1000:"
                  " at least 10^28 consecutive certified binary orders.")



if __name__ == '__main__':
    verify()
