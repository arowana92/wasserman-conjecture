"""Independent forward binary and factorial-floor checks for five related weights."""
from argparse import ArgumentParser
from bisect import bisect_left
from collections import deque
from functools import reduce
from math import comb, lcm
from operator import or_


def cut(w, a):
    # Each tuple represents (slope*t + offset)/denominator.
    return [(a, 0, 1), (w-a, 1, 1)]


def binary_residues(rows, low, width, modulus):
    forms = [form for row in rows for form in row]
    ends = [0]
    for row in rows:
        ends.append(ends[-1]+len(row))

    def advance(state, bit):
        numerators = [c+a*bit for c, (a, _, _) in zip(state, forms)]
        digits = [n % 2 for n in numerators]  # Every denominator is odd.
        if any(sum(digits[i:j]) > 1 for i, j in zip(ends, ends[1:])):
            return None
        return tuple((n-d*e)//2 for n, e, (_, _, d) in zip(numerators, digits, forms))

    def accepts(state):
        high = []
        for c, (_, _, d) in zip(state, forms):
            if c < 0 or c % d:
                return False
            high.append(c//d)
        return all(sum(high[i:j]) == reduce(or_, high[i:j])
                   for i, j in zip(ends, ends[1:]))

    state = tuple(b for _, b, _ in forms)
    for i in range(width):
        state = advance(state, (low >> i) & 1)
        if state is None:
            return set()
    initial = state, low % modulus, 2**width % modulus
    todo = deque([initial])
    seen = {initial}
    residues = set()
    while todo:
        state, value, place = todo.popleft()
        if accepts(state):
            residues.add(value)
        for bit in (0, 1):
            child = advance(state, bit)
            key = child, (value+bit*place) % modulus, 2*place % modulus
            if child is not None and key not in seen:
                seen.add(key)
                todo.append(key)
    return residues


def assignments():
    found = []
    for w in (8, 16):
        for a in range(3, w, 2):
            if w == 8 and a == 7:
                continue
            low = pow(a, -1, w)
            for u in range(3, a+1, 2):
                if a % u:
                    continue
                for v, b in ((7, 6), (9, 4), (9, 7), (9, 8)):
                    if (v == 7 or b == 4) and low % 8 != 3:
                        continue
                    if b == 8 and low % 4 != 3:
                        continue
                    modulus = lcm(u, v, 9 if v == 7 else 1)
                    other = ([(9*y, y, 7) for y in (1, 2, 4)]
                             if (v, b) == (9, 7) else cut(v, b))
                    seen = binary_residues([cut(w, a), other], low, w.bit_length()-1, modulus)
                    if any((w*t+1) % u == 0 and (b*t-1) % v == 0
                           and (v != 7 or t % 9 == 2) for t in seen):
                        found.append((w, a, u, v, b))
    assert found == [(16, 3, 3, 7, 6), (16, 11, 11, 9, 4)]
    for a in (3, 5):
        b = a+8
        seen = binary_residues([cut(8, a), cut(16, b)], pow(b, -1, 16), 4, a)
        assert all((8*t+1) % a for t in seen)
    print("Forward rational binary automata confirm all assignment exclusions")


def floor_search(cuts, limit, modulus, residue):
    layers = []
    for p, w, a in cuts:
        powers = []
        q = p
        while q <= w*limit+1:
            powers.append(q)
            q *= p
        layers.append((powers, w, a))
    stack = [(2, limit)]
    visited = 0
    hits = []
    while stack:
        low, high = stack.pop()
        low += (residue-low) % modulus
        high -= (high-residue) % modulus
        if low > high:
            continue
        visited += 1
        excluded = False
        for powers, w, a in layers:
            n, x, y = w*low+1, a*high, (w-a)*high+1
            # Skipping small layers only weakens a rejection. At a singleton
            # every layer is considered, giving the exact binomial valuation.
            start = bisect_left(powers, high-low)
            for i in range(start, len(powers)):
                q = powers[i]
                if n//q > x//q+y//q:
                    excluded = True
                    break
            if excluded:
                break
        if excluded:
            continue
        if low == high:
            hits.append(low)
        else:
            mid = (low+high)//2
            stack.extend(((mid+1, high), (low, mid)))
    return visited, hits


def finite_check(digits):
    for p, w, a in ((2, 8, 3), (3, 9, 4), (7, 7, 6), (11, 16, 11)):
        actual = floor_search([(p, w, a)], 140, 1, 0)[1]
        assert actual == [t for t in range(2, 141) if comb(w*t+1, a*t) % p]
        assert actual
    for v, cuts, modulus, residue in (
        (7, [(2, 16, 3), (2, 7, 6), (3, 7, 6), (3, 16, 3), (7, 7, 6)], 1008, 587),
        (9, [(2, 16, 11), (2, 9, 4), (3, 9, 4), (11, 16, 11)], 1584, 1267),
    ):
        if v == 7:
            assert residue % 16 == residue % 72 == 11 and 6*residue % 7 == 1
        else:
            assert residue % 144 == 115 and (16*residue+1) % 11 == 0
        visited, hits = floor_search(cuts, 10**digits, modulus, residue)
        print(f"Pair 16 and {v}: {visited:,} nonempty intervals, {len(hits)} survivors "
              f"through 10^{digits}", flush=True)
        assert not hits


if __name__ == "__main__":
    cli = ArgumentParser(description=__doc__)
    cli.add_argument("--digits", type=int, default=1000)
    args = cli.parse_args()
    assignments()
    finite_check(args.digits)
