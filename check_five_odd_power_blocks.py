"""Cross-check both binary traversals and finish their arbitrary-precision tails."""

import argparse
from pathlib import Path
import subprocess


def read_run(program, extra):
    result = subprocess.check_output([program, "1000", *extra], text=True)
    levels = {}
    frontier = {}
    for line in result.splitlines():
        label, *fields = line.split()
        key, *values = map(int, fields)
        if label == "level":
            assert key not in levels and len(values) == 1
            levels[key] = values[0]
        else:
            assert label == "suffix" and key not in frontier
            frontier[key] = values
    return levels, frontier


def disjoint(value, exponent, bits):
    modulus = 1 << bits
    q = pow(value, exponent, modulus)
    return (3*q & 2*q & (modulus-1)) == 0


def finish(program):
    levels, frontier = read_run(program, [])
    assert (levels, frontier) == read_run(program, ["--logs"])
    assert sum(levels.values()) == 73_328_239 and levels[120] == len(frontier) == 18
    assert max(levels.items(), key=lambda item:item[1]) == (30,4_918_691)
    for r, exponents in frontier.items():
        assert exponents == [e for e in range(1,2000,2) if disjoint(r,e,120)]
    bits = 120
    while frontier:
        new = {}
        for r, exponents in frontier.items():
            for x in (r, r+(1 << bits)):
                retained = [e for e in exponents if disjoint(x,e,bits+1)]
                if len(retained) >= 3:
                    new[x] = retained
        bits += 1
        frontier = new
    assert bits == 139
    print("Both coordinate systems agree at every level through 120 bits.")
    print("All 18 remaining prefixes terminate by 139 bits.")


def small_certificates():
    for a,b,precision in [(1,19,9),(1,121,15),(1,251,14),(121,251,9)]:
        def compatible(bits):
            return any(disjoint(q,a,bits) and disjoint(q,b,bits)
                       for q in range(3,1 << bits,8))
        assert compatible(precision-1) and not compatible(precision)
    m = 128
    s = [e for e in range(1,m,2) if disjoint(3,e,9)]
    ratios = {a*pow(b,-1,m) % m for a in s for b in s}
    assert len(s) == 16 and set(range(1,m,2))-ratios == {19,27}
    print("Reproduced the nine-bit obstruction and all three pair exclusions.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    args = parser.parse_args()
    small_certificates()
    finish(args.executable.resolve())
