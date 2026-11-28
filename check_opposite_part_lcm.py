"""Test the four-factor inequality on pairs passing all full endpoint powers."""
import argparse
from functools import lru_cache
from itertools import product
from math import lcm, prod


@lru_cache(None)
def blocks(n):
    result = []
    p = 2
    while p*p <= n:
        power = 1
        while n%p == 0:
            n //= p
            power *= p
        if power > 1:
            result.append(power)
        p = 3 if p == 2 else p+2
    if n > 1:
        result.append(n)
    return tuple(result)


def minimum_cost(m):
    groups = {(1, 1, 1, 1)}
    for power in blocks(m):
        groups = {tuple(sorted(group[:i]+(group[i]*power,)+group[i+1:]))
                  for group in groups for i in range(4)}
    # Pairing the two smallest factors together minimizes the product of sums.
    return min((a+b)*(c+d) for a,b,c,d in groups)


def distribute(values):
    for bits in range(1 << len(values)):
        yield (tuple(v for i,v in enumerate(values) if bits >> i & 1),
               tuple(v for i,v in enumerate(values) if not (bits >> i & 1)))


def sample_rows(n, scale, assigned):
    size, modulus = n//scale, prod(assigned)
    steps = [modulus//q * pow(scale*(modulus//q), -1, q) for q in assigned]
    rows = set()
    for positions in product(range(3), repeat=len(steps)):
        least = [0, 0, 0]
        for step, i in zip(steps, positions):
            least[i] += step
        least = [((a-1)%modulus)+1 for a in least]
        spare = size-sum(least)
        if spare < 0:
            continue
        assert spare%modulus == 0
        # Choose one representative per ordered residue pattern.
        least[-1] += spare
        rows.add(tuple(sorted(scale*a for a in least)))
    return rows


def check_pair(n, first, second):
    assert sum(first) == sum(second) == n
    for endpoint in (n, n-1):
        for q in blocks(endpoint):
            assert min(sum(a%q for a in row) for row in (first,second)) == n%q
    minimum = n*n
    for a in first:
        for b in second:
            m, remainder = divmod(lcm(a,b), n)
            assert not remainder
            k, remainder = divmod(m*(n-a)*(n-b), n*(n-1))
            assert not remainder and minimum_cost(m) <= k < m
            minimum = min(minimum, m)
    return minimum


def scan(limit):
    count = 0
    weights = set()
    smallest = None
    for n in range(3, limit+1):
        powers = blocks(n-1)
        if len(powers) < 3:
            continue
        pairs = set()
        for first, second in distribute(blocks(n)):
            u, v = prod(first), prod(second)
            if not 3 <= u <= v:
                continue
            for s,t in distribute(powers):
                if any(q >= v for q in s) or any(q >= u for q in t):
                    # At least two normalized parts are positive multiples of q.
                    continue
                candidates = sample_rows(n, u, s)
                if candidates:
                    pairs.update(product(candidates, sample_rows(n, v, t)))
        for first, second in pairs:
            value = check_pair(n, first, second)
            item = (value, n, first, second)
            smallest = item if smallest is None else min(smallest, item)
            count += 1
            weights.add(n)
    print("Distinct weights:", len(weights), "row pairs:", count,
          "part pairs:", 9*count, "least quotient and example:", smallest)
    if limit == 10000:
        assert (len(weights), count, smallest) == (35, 41, (35, 3256, (280, 1240, 1736),
                                                           (814, 1221, 1221)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("limit", type=int, nargs="?", default=10000)
    options = parser.parse_args()
    scan(options.limit)
