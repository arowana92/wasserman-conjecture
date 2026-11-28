"""Exact unrestricted and endpoint-only certificate sizes at weight 12,189,629."""

from itertools import combinations
from math import isqrt, prod
from pathlib import Path

from search_endpoint import find_row


N = 12189629
WEIGHT = (17, 53, 83, 163)
PREDECESSOR = (2, 11, 29, 41, 233)
ENDPOINT = tuple(sorted(WEIGHT + PREDECESSOR))
SEVEN = (17, 29, 41, 53, 83, 163, 233)
NEAR = 12189623


def valuation(n, k, p):
    """Legendre valuation of the binomial n choose k."""
    value = 0
    w = p
    while w <= n:
        value += n // w - k // w - (n - k) // w
        w *= p
    return value


def row_valuation(n, a, b, p):
    return valuation(n, a, p) + valuation(n - a, b, p)


def lower_bound():
    path = Path(__file__).with_name("six_prime_witnesses.txt")
    supports = []
    labels = set()
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        p, q, a, b = map(int, line.split())
        if (p, q) in labels or not (p < q and min(a, b, N - a - b) > 0):
            raise ArithmeticError(line)
        labels.add((p, q))
        avoided = {s for s in ENDPOINT if row_valuation(N, a, b, s) == 0}
        if not {p, q} <= avoided:
            raise ArithmeticError(("invalid witness", line))
        supports.append(avoided)
    missing = set(combinations(ENDPOINT, 2)) - labels
    if missing != {(17, 163), (83, 163)}:
        raise ArithmeticError(("unexpected pairs without witnesses", missing))
    pairs = [a | b for a in supports for b in supports]
    for group in combinations(ENDPOINT, 4):
        if not any(set(group) <= union for union in pairs):
            raise ArithmeticError(group)
    return len(supports)


def common_cuts(group):
    # Use the first residue layer at each endpoint prime.
    residues = [0]
    period = 1
    for p in group:
        target = (0,) if p in WEIGHT else (0, 1)
        residues = [
            r + period * ((v - r) * pow(period, -1, p) % p)
            for r in residues
            for v in target
        ]
        period *= p
    return sorted(
        a
        for r in residues
        for a in range(r or period, N - 1, period)
        if all(valuation(N, a, p) == 0 for p in group)
    )


def upper_bound():
    checks = 0
    for group in combinations(SEVEN, 3):
        cuts = common_cuts(group)
        for i, a in enumerate(cuts):
            if 3 * a > N:
                break
            for b in cuts[i:]:
                if a + 2 * b > N:
                    break
                checks += 1
                if all(valuation(N - a, b, p) == 0 for p in group):
                    raise ArithmeticError(("simultaneously avoided", group, a, b))
        # Independent CRT on full rows, with direct base-p digit additions.
        wa = [(p, p) for p in group if p in WEIGHT]
        ma = [(p, p) for p in group if p in PREDECESSOR]
        if (
            find_row(N, prod(p for p, _ in wa), wa, prod(p for p, _ in ma), ma)
            is not None
        ):
            raise ArithmeticError(("CRT mismatch", group))
    return checks


if __name__ == "__main__":
    if any(any(p % d == 0 for d in range(2, isqrt(p) + 1)) for p in ENDPOINT + (NEAR,)):
        raise ArithmeticError("composite base")
    if prod(WEIGHT) != N or 2 * prod(PREDECESSOR) != N - 1 or N - NEAR != 6:
        raise ArithmeticError("factorization or prime gap")
    print(
        f"Verified {lower_bound()} pair witnesses and all 126 endpoint subsets of size four."
    )
    print(
        f"Excluded all 35 avoided triples of the seven primes; {upper_bound()} binomial comparisons."
    )
    near_rows = [(a, b) for a in range(1, 7) for b in range(a, 7 - a)]
    if not all(row_valuation(N, a, b, p) > 0 for a, b in near_rows for p in SEVEN):
        raise ArithmeticError("near-prime row")
    print(
        f"All {len(near_rows)} rows avoiding {NEAR} carry at each of the seven endpoint primes."
    )
    print(
        "Minimum sizes: endpoint pairs 5, arbitrary triples 6, endpoint-only triples 7."
    )
