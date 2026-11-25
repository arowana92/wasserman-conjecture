"""Generate the short prime lists and check their performance through 5000."""

from argparse import ArgumentParser
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

from flint import fmpz
from prime_certificates import find_cover, free_count
from search_even_endpoint import prime_powers


@lru_cache(maxsize=8)
def retained(n):
    if n <= 1:
        return ()
    small = 1
    selected = []
    for p, w in sorted(prime_powers(n).items()):
        small *= w
        if small**25 > n:
            selected.append(p)
    return tuple(selected)


def certificate_list(n):
    if n < 3:
        raise ValueError("weight must be at least three")
    r = n if n % 2 else n - 1
    while not fmpz(r).is_prime():
        r -= 2
    primes = {r} | {p for p in (2, 3, 5) if p <= n}
    for h in range(6):
        primes.update(retained(n-h))
    return sorted(primes)


def check_digit_bounds():
    for p in (2, 3, 5):
        f = sum(Fraction(1, 3)**d for d in range(p))
        g = sum(Fraction(1, 2)**d for d in range(p))
        exponent = 50*(p-1) + 1
        if f**300 * 3**exponent >= p**213 or g**200 * 2**exponent >= p**172:
            raise ArithmeticError(("digit bound", p))


def check_table():
    # The existing table was verified by exact factorial valuations on all rows.
    path = Path(__file__).with_name("certificates_5000.txt")
    covered = size_sum = maximum = additional = 0
    failures = []
    for line in path.read_text().splitlines():
        n, *known = map(int, line.split())
        chosen = certificate_list(n)
        size_sum += len(chosen)
        maximum = max(maximum, len(chosen))
        if not set(known).issubset(chosen):
            ordered = sorted(chosen, key=lambda p: free_count(n, p))
            if find_cover(n, ordered) is not None:
                failures.append(n)
                continue
            additional += 1
        covered += 1
    print(f"Certified {covered} weights, {additional} by a separate prime-set search.")
    print(f"Total listed primes {size_sum}, maximum list size {maximum}.")
    print(f"Weights at which this list fails: {failures}")


if __name__ == "__main__":
    cli = ArgumentParser(description=__doc__)
    cli.add_argument("weight", nargs="?", type=int, help="print the candidate prime list")
    args = cli.parse_args()
    check_digit_bounds()
    if args.weight is None:
        check_table()
    else:
        print(*certificate_list(args.weight))
