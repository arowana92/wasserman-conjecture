"""Check shift-independent exclusions at 1+3^e (3^a+2)(3^a+4)."""

from argparse import ArgumentParser
from collections import Counter
from functools import lru_cache
from itertools import combinations
from math import gcd, prod
from pathlib import Path

from flint import fmpz

from check_six_bit_rays import normalized_rows
from prime_certificates import primes_through
from search_even_endpoint import prime_powers


CERTIFICATE = Path(__file__).with_name("ternary_quadratic_factors.txt")


def full(n, q):
    if q < 2 or n % q or not fmpz(q).is_prime():
        raise ArithmeticError(("invalid prime divisor", n, q))
    w = q
    while n % (w * q) == 0:
        w *= q
    return w


def first_carries(s, q, w):
    t = s // w
    level = q
    while level <= w:
        if (-t) % level + (-t * pow(3, -1, level)) % level >= level:
            return True
        level *= q
    return False


def required_companion_prime(q):
    # 2 must have odd multiplicative order at every prime taken by A here.
    odd = (q - 1) // ((q - 1) & -(q - 1))
    return pow(2, odd, q) != 1


def check(a, lists, counts):
    X = 3**a
    s, D, M = X + 1, X + 2, X + 4
    bad_u, bad_b, crt = lists
    if any(len(set(group)) != len(group) for group in lists):
        raise ArithmeticError("repeated factor")
    upper, lower = s, M
    for q in bad_u:
        w = full(s, q)
        if not first_carries(s, q, w):
            raise ArithmeticError(("missing weight carry", a, q))
        upper //= w
    for q in bad_b:
        w = full(D, q)
        if not required_companion_prime(q):
            raise ArithmeticError(("odd order", a, q))
        lower *= w
    if 27 * lower**2 > upper**3:
        counts["size"] += 1
        return True
    moduli = {q: full(M, q) for q in crt}
    if not moduli or prod(moduli.values()) != M:
        return False
    if any(2 * w + 1 > upper for w in moduli.values()):
        counts["single power"] += 1
        return True
    for z in normalized_rows(s, tuple(moduli.values())):
        counts["CRT rows"] += 1
        if upper % (s // gcd(*z)):
            continue
        if any(gcd((x + y) // gcd(x, y), M) > 1 for x, y in combinations(z, 2)):
            continue
        if any(
            sum((-(x // w) * pow(3, -1, w)) % w for x in z if x % w == 0) >= w
            for w in moduli.values()
        ):
            continue
        return False
    counts["complete CRT"] += 1
    return True


def read_certificate(path=CERTIFICATE):
    records = {}
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        index, *fields = line.split()
        lists = tuple([] if f == "-" else list(map(int, f.split(","))) for f in fields)
        if len(lists) != 3 or int(index) in records:
            raise ArithmeticError(("invalid record", line))
        records[int(index)] = lists
    return records


@lru_cache(maxsize=2048)
def split_small(n):
    found = []
    for q in primes_through(10000):
        if n % q == 0:
            found.append(q)
            while n % q == 0:
                n //= q
    if n > 1 and fmpz(n).is_prime():
        found.append(n)
        n = 1
    return found, n


@lru_cache(maxsize=1024)
def cyclotomic_at_three(n):
    from sympy import divisors

    return (3**n - 1) // prod(cyclotomic_at_three(d) for d in divisors(n) if d < n)


def make(a, deep):
    from sympy import divisors

    X = 3**a
    s, D, M = X + 1, X + 2, X + 4
    components = [cyclotomic_at_three(d) for d in divisors(2 * a) if a % d]
    pieces = [M]
    if a % 4 == 0:
        y = 3 ** (a // 4)
        pieces = [y * y - 2 * y + 2, y * y + 2 * y + 2]
    todo = []
    groups = [set(), set(), set()]
    for group, values in enumerate((components, [D], pieces)):
        for n in values:
            small, rest = split_small(n)
            groups[group].update(small)
            if rest > 1:
                todo.append((rest, group))

    def attempt():
        bad_u = [q for q in groups[0] if first_carries(s, q, full(s, q))]
        bad_b = [q for q in groups[1] if required_companion_prime(q)]
        crt = sorted(groups[2])
        if prod(full(M, q) for q in crt) != M:
            crt = []
        lists = [sorted(bad_u), sorted(bad_b), crt]
        if not check(a, lists, Counter()):
            return None
        # Retain a small certificate; preserve complete CRT factorizations.
        trial = [lists[0][:], lists[1][:], []]
        if check(a, trial, Counter()):
            lists = trial
        for i in (0, 1):
            for q in lists[i][:]:
                lists[i].remove(q)
                if not check(a, lists, Counter()):
                    lists[i].append(q)
            lists[i].sort()
        return lists

    result = attempt()
    if not deep or result is not None:
        return result
    for n, group in sorted(todo):
        groups[group].update(prime_powers(n))
        result = attempt()
        if result is not None:
            return result
    return None


def encode(a, groups):
    return " ".join([str(a), *(",".join(map(str, g)) or "-" for g in groups)])


def check_outer_limit(limit, start):
    primes = primes_through((limit + 1) // start)
    a = start
    while True:
        bound = (limit + 1) // a
        F = 3**a + 2
        # The left side decreases with a and the right side increases.
        if 1 << (2 * bound) <= 27 * F:
            break
        powers = []
        for q in primes:
            if q > bound:
                break
            if q == 2:
                continue
            w = q
            while w * q <= bound:
                w *= q
            powers.append(w)
        multiple = fmpz(prod(powers))
        cap = (pow(fmpz(2), multiple, fmpz(F)) - 1).gcd(F)
        if cap**2 > 27 * F:
            raise ArithmeticError(("inconclusive aggregate gcd", a, cap))
        a += 1
    print(
        f"Every outer exponent through {limit} certified; {a-start} aggregate gcds from a={start}."
    )


if __name__ == "__main__":
    cli = ArgumentParser(description=__doc__)
    cli.add_argument("--generate", type=int, metavar="LIMIT")
    cli.add_argument("--deep", action="store_true")
    cli.add_argument("--start", type=int, default=1)
    cli.add_argument("--outer-limit", type=int)
    args = cli.parse_args()
    if args.generate:
        for a in range(args.start, args.generate + 1):
            groups = make(a, args.deep)
            print(
                encode(a, groups) if groups is not None else f"# unresolved {a}",
                flush=True,
            )
    else:
        records = read_certificate()
        if sorted(records) != list(range(1, max(records) + 1)):
            raise ArithmeticError("nonconsecutive certificate")
        counts = Counter()
        for a, lists in records.items():
            if not check(a, lists, counts):
                raise ArithmeticError(("inconclusive", a))
        print(
            f"All outer exponents certified for inner exponents 1 through {max(records)}: {dict(counts)}"
        )
        if args.outer_limit is not None:
            check_outer_limit(args.outer_limit, max(records) + 1)
