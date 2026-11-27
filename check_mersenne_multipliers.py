"""Check all binary outer exponents simultaneously at all Mersenne primes."""

from argparse import ArgumentParser
from collections import Counter, deque
from functools import lru_cache
from math import gcd

from search_endpoint import atoms, smallest_factors


@lru_cache(None)
def subgroup(modulus):
    value = 1 % modulus
    powers = set()
    while value not in powers:
        powers.add(value)
        value = 2 * value % modulus
    return powers


def shapes(d, powers):
    for i in range(1, d // 4 + 1):
        for j in range(i + 1, (d - 1) // 2 - i + 1):
            if (i & -i) == (j & -j):
                continue
            bound = (d - 1) // (i + j)
            for a in range(1, bound // 2 + 1):
                for b in range(a, bound - a + 1):
                    if all(any(z % w == 0 for z in (i, j, a, b)) for w in powers):
                        yield i, j, a, b


def assignments(powers, i, j, a, b):
    choices = {(1, 1)}
    for w in powers:
        if (i % w == j % w == 0) or (a % w == b % w == 0):
            continue
        next_choices = set()
        for first, second in choices:
            if (i % w == 0 or j % w == 0) and gcd(i + j, w) == 1:
                next_choices.add((first * w, second))
            if (a % w == 0 or b % w == 0) and gcd(a + b, w) == 1:
                next_choices.add((first, second * w))
        choices = next_choices
    return choices


def admissible(d, shape, assigned, fixed_height=None):
    i, j, a, b = shape
    first, second = assigned
    modulus = first * second
    half = pow(2, -1, modulus)
    offset = (a + b) * first * pow(first, -1, second) % modulus if second > 1 else 0
    allowed = subgroup(first)
    length = (2 * d - 1).bit_length()

    def next_state(state, digit, target):
        x, y, residue = state
        x += i * digit
        y += j * digit
        if (x & 1) + (y & 1) <= target:
            return x // 2, y // 2, (residue + digit) * half % modulus
        return None

    @lru_cache(None)
    def high_inputs(x, y):
        values = []
        for v in range(a + b, (d - 1 - x - y) // (i + j) + 1):
            r, s = i * v + x, j * v + y
            if not ((r & s) or ((r | s) & ~(d - 1))):
                values.append(v)
        return values

    # For a fraction t in [0,1), (x,y)=(floor(i*t),floor(j*t)).
    # As i<j, each interval for y meets at most two intervals for x.
    if not any(high_inputs(x, y)
               for y in range(j) for x in {i * y // j, (i * (y + 1) - 1) // j}):
        return False

    def accepts(state):
        x, y, residue = state
        return any((v + residue - a - b) % second == 0
                   and (i + j) * (v + residue) % first in allowed for v in high_inputs(x, y))

    def advance(states, digits, target):
        return {new for state in states for digit in digits
                if (new := next_state(state, digit, target)) is not None}

    states = {(0, 0, offset)}
    for k in range(length):
        states = advance(states, (1,) if k == 0 else (0, 1), (-d >> k) & 1)
        if not states:
            return False

    if fixed_height is not None:
        assert fixed_height >= length
        for _ in range(length, fixed_height):
            states = advance(states, (0, 1), 1)
        return any(map(accepts, states))

    pending = deque(states)
    while pending:
        state = pending.popleft()
        if accepts(state):
            return True
        for new in advance((state,), (0, 1), 1):
            if new not in states:
                states.add(new)
                pending.append(new)
    return False


def direct_at_height(d, shape, assigned, h):
    i, j, a, b = shape
    first, second = assigned
    r = 2**h - 1
    total = d * r
    allowed = subgroup(first)
    for u in range((r * (a + b) + 1) | 1, total // (i + j) + 1, 2):
        x, y = i * u, j * u
        if (x & y) or ((x | y) & ~total):
            continue
        if (u - r * (a + b)) % second == 0 and u * (i + j) % first in allowed:
            return True
    return False


def main(limit, cross_limit):
    spf = smallest_factors(limit)
    residual = set()
    counts = Counter()
    for d in range(1, limit + 1, 2):
        powers = [w for _, w in atoms(d, spf)]
        for shape in shapes(d, powers):
            counts["shapes"] += 1
            for assigned in assignments(powers, *shape):
                counts["assignments"] += 1
                works = admissible(d, shape, assigned)
                if works:
                    residual.add(d)
                    counts["surviving assignments"] += 1
                if d <= cross_limit:
                    start = (2 * d - 1).bit_length()
                    for h in range(start, start + 3):
                        finite = admissible(d, shape, assigned, fixed_height=h)
                        assert finite == direct_at_height(d, shape, assigned, h)
                        assert not finite or works
                        counts["direct fixed-height comparisons"] += 1
    good = sum(d // (d & -d) not in residual for d in range(1, limit + 1))
    print(dict(counts))
    print("Odd multipliers still allowed:", sorted(residual))
    print("Certified multipliers including even ones:", good, "of", limit)
    return residual


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--limit", type=int, default=512)
    parser.add_argument("--cross-limit", type=int, default=127)
    args = parser.parse_args()
    main(args.limit, args.cross_limit)
