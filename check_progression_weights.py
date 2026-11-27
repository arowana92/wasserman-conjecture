"""Exact controls for arithmetic progression multiplicities and entropy constants."""
from collections import Counter
from fractions import Fraction as F
from math import gcd

from check_binary_windows import LN2_HI, LN2_LO, ln_bounds


def entropy_bounds(x):
    lo = hi = 0
    for mass in (x, 1-x):
        l, r = ln_bounds(1/mass)
        lo += mass*l/LN2_HI
        hi += mass*r/LN2_LO
    return lo, hi


def check_constants():
    alo, ahi = entropy_bounds(F(1, 6))
    blo, bhi = entropy_bounds(F(1, 3))
    assert F(65002242, 10**8) < alo <= ahi < F(65002243, 10**8)
    assert F(72405194, 10**8) < blo/(1+blo-alo)
    assert bhi/(1+bhi-ahi) < F(72405195, 10**8)
    assert alo*(1-ahi) > F(2, 9) and ahi < F(2, 3) < blo


def check_blocks():
    runs = shared = 0
    for step in range(1, 65):
        u = (step-1).bit_length()
        t = (step & -step).bit_length()-1
        for length in (65, 96, 127, 128, 189):
            l = length.bit_length()-1
            for start in (256, 257, 511, 512, 777):
                weights = range(start, start+length*step, step)
                low = Counter((n >> t) % (1 << (u-t)) for n in weights)
                high = Counter((n >> u) % (1 << l) for n in weights)
                assert max(low.values()) <= (length+(1 << (u-t))-1) >> (u-t)
                assert max(high.values()) <= 4
                assert len({n >> (l+u) for n in weights}) <= 3
                for modulus in (12, 30, 64, 105, 231):
                    f = gcd(step, modulus)
                    groups = Counter(n % modulus for n in weights)
                    assert max(groups.values()) <= (length*f+modulus-1)//modulus
                    shared += f > 1
                runs += 1
    for f in range(1, 129):
        bits = f.bit_length()-1
        keys = {(z >> bits, z % f) for z in range(1024)}
        assert len(keys) == 1024
    print("Progressions:", runs, "noncoprime modulus tests:", shared)


if __name__ == "__main__":
    check_constants()
    check_blocks()
    print("Both threshold enclosures and all block tests passed.")
