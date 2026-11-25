"""Find and verify optimal four-prime certificates at adjacent weights."""

from argparse import ArgumentParser
from bisect import bisect_right
from collections import Counter
from math import comb, isqrt

from check_minimality import digit_sum
from prime_certificates import primes_through


def prime_divisors(n, trial):
    divisors = []
    for p in trial:
        if p*p > n:
            break
        if n % p == 0:
            divisors.append(p)
            while n % p == 0:
                n //= p
    if n > 1:
        divisors.append(n)
    return divisors


def previous_prime(n, trial):
    candidate = n if n % 2 else n-1
    while candidate > 2:
        if all(candidate % p for p in trial[:bisect_right(trial, isqrt(candidate))]):
            return candidate
        candidate -= 2
    return 2


def inspect_pair(p, q, trial, cutoff):
    n = p*q
    before = prime_divisors(n-1, trial)
    after = prime_divisors(n+1, trial)
    minus, plus = before[-1], after[-1]
    if min(minus, plus) <= cutoff:
        return minus, plus, None
    r = previous_prime(n, trial)
    # Use the same near prime for both weights and the stronger common gap bound.
    if 2*r <= n+1 or n+1-r >= min(p, q, minus, plus):
        return minus, plus, None
    if any(digit_sum(n+1, base) <= 2 for base in [p, q]+after):
        return minus, plus, None
    return minus, plus, ((p, q, minus, r), (p, q, plus, r))


def scan(scale):
    trial = primes_through(2*scale+2)
    interval = [p for p in trial if scale < p and 3*p <= 4*scale]
    counts = Counter()
    first = None
    first_primes = [p for p in interval if p % 4 == 1]
    second_primes = [q for q in interval if q % 4 == 3]
    for p in first_primes:
        for q in second_primes:
            minus, plus, certs = inspect_pair(p, q, trial, 2*scale)
            counts["pairs"] += 1
            counts["large minus"] += minus > 2*scale
            counts["large plus"] += plus > 2*scale
            counts["both large"] += min(minus, plus) > 2*scale
            if certs:
                counts["optimal adjacent pairs"] += 1
                if first is None:
                    first = p*q, certs
    print(scale, dict(counts), "first", first)
    return counts


def exact(n, certificate):
    if len(set(certificate)) != 4 or any(
        p < 2 or any(p % d == 0 for d in range(2, isqrt(p)+1)) for p in certificate
    ):
        raise ArithmeticError(("invalid prime list", certificate))
    counts = Counter()
    for a in range(1, n//3+1):
        first = comb(n, a)
        for b in range(a, (n-a)//2+1):
            value = first*comb(n-a, b)
            missed = tuple(p for p in certificate if value % p)
            if len(missed) > 1:
                raise ArithmeticError((n, a, b, missed))
            counts[missed] += 1
    # Test every possible universal prime directly on digit sums.
    if any(digit_sum(n, p) <= 2 for p in primes_through(n)):
        raise ArithmeticError(("universal prime", n))
    print(f"Weight {n}, primes {certificate}, {sum(counts.values())} exact rows, "
          f"rows avoiding each prime {[counts[p,] for p in certificate]}.")


if __name__ == "__main__":
    cli = ArgumentParser(description=__doc__)
    cli.add_argument("scales", type=int, nargs="*", default=[100, 1000, 10000])
    args = cli.parse_args()
    _, _, example = inspect_pair(17, 11, primes_through(30), 17)
    if example is None:
        raise ArithmeticError("Missing adjacent certificate")
    for n, certificate in zip((187, 188), example):
        exact(n, certificate)
    for scale in args.scales:
        if scale < 3:
            cli.error("use scales of at least 3")
        scan(scale)
