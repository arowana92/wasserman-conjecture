"""List binary offsets passing the necessary test at unbounded heights."""

from argparse import ArgumentParser
from functools import lru_cache
from itertools import combinations_with_replacement
from math import gcd

from check_binary_gap_bounds import entropy_below


def mark_offsets(limit):
    possible = bytearray(limit+1)
    for t in [1, *range(3, limit//3+1, 2)]:
        values = [-1] if t == 1 else []
        odd = t
        while odd <= limit:
            z = odd
            while z <= limit:
                values.append(z)
                z *= 2
            if t == 1:
                break
            odd *= t
        values.sort()
        for x, y, z in combinations_with_replacement(values, 3):
            total = x+y+z
            if total > limit or y < 0:
                continue
            if x >= 0 and (x | y | z) == total:
                possible[total] = 1
            if total < limit and (
                ((x+1) | y | z) == total+1 or
                (x >= 0 and ((x | (y+1) | z) == total+1 or (x | y | (z+1)) == total+1))
            ):
                possible[total+1] = 1
    return possible


@lru_cache(maxsize=None)
def core(z):
    z = abs(z)
    z //= z & -z
    if z == 1:
        return 1
    factors = []
    p = 3
    while p*p <= z:
        e = 0
        while z % p == 0:
            z //= p
            e += 1
        if e:
            factors.append((p,e))
        p += 2
    if z > 1:
        factors.append((z,1))
    power = gcd(*(e for _,e in factors))
    ans = 1
    for p,e in factors:
        ans *= p**(e//power)
    return ans


def by_digits(c):
    def subsets(x):
        y = x
        while True:
            yield y
            if y == 0:
                return
            y = (y-1) & x
    for a in subsets(c):
        remaining = c ^ a
        for b in subsets(remaining):
            row = (a, b, remaining ^ b)
            for i in (-1,0,1,2):
                z = [v-int(i==j) for j,v in enumerate(row)]
                if all(z) and len({core(v) for v in z}) == 1:
                    return True
    return False


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--limit", type=int, default=2**20)
    parser.add_argument("--cross", type=int, default=4096)
    args = parser.parse_args()
    assert entropy_below(7,10,89,100)
    assert entropy_below(7,30,79,100)
    assert entropy_below(1,100,81,1000)
    possible = mark_offsets(args.limit)
    for c in range(2, min(args.limit,args.cross)+1):
        assert bool(possible[c]) == by_digits(c), c
    for cutoff in sorted({v for v in (256,4096,65536,args.limit) if v<=args.limit}):
        passing = sum(possible[2:cutoff+1])
        print(f"2..{cutoff}: {passing} pass the necessary test, {cutoff-1-passing} excluded")


if __name__ == "__main__":
    main()
