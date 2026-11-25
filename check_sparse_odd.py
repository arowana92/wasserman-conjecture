"""Verify odd four- and five-bit weights by digit partitions and factorial bounds."""

from argparse import ArgumentParser
from collections import Counter
from itertools import combinations
from math import comb, gcd, prod

from flint import fmpz

from check_three_bit_range import exclude_without_endpoint


def partitions(units):
    """Generate unordered partitions of the digit units into at most three parts."""
    if not units:
        yield ()
        return
    digit = units[-1]
    for row in partitions(units[:-1]):
        if len(row) < 3:
            yield (*row, digit)
        for i, part in enumerate(row):
            yield row[:i] + (part+digit,) + row[i+1:]


def gap_suffices(n):
    r = n
    while not fmpz(r).is_prime():
        r -= 2
    g = n - r
    if g <= 10:
        return True
    lower = prod(((n-h)//gcd(n-h, 60))**(6-h) for h in range(6))
    return 3**18 * lower > 64*n**18 * g**15


def check(bits, digits):
    counts = Counter()
    for positions in combinations(range(1, bits), digits-1):
        units = (1,) + tuple(1 << i for i in positions)
        n = sum(units)
        if n <= 1_000_000:
            counts["previous million-weight theorem"] += 1
        elif gap_suffices(n):
            counts["gap and height six"] += 1
        else:
            counts["weights requiring digit partitions"] += 1
            rows = [sorted(t) for t in partitions(units) if len(t) == 3]
            assert len(rows) == (3**digits - 3*2**digits + 3)//6
            for row in rows:
                exclude_without_endpoint(n, row, counts)
    print(f"Odd weights below 2^{bits} with {digits} binary ones: {comb(bits-1, digits-1)}")
    print(dict(counts), flush=True)


if __name__ == "__main__":
    cli = ArgumentParser(description=__doc__)
    cli.add_argument("--four-bits", type=int, default=128, metavar="POSITIONS")
    cli.add_argument("--five-bits", type=int, default=64, metavar="POSITIONS")
    args = cli.parse_args()
    check(args.four_bits, 4)
    check(args.five_bits, 5)
