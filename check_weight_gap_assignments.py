"""Independently test whether two normalized triples can cover all residue constraints."""

from argparse import ArgumentParser
from functools import lru_cache
from itertools import product
from math import prod

from check_predecessor_weight_gap import criterion, factorizations
from check_three_factors import valuation


def permits_row(total, moduli, scale=1, bases=()):
    if total < 3 or any(2*w + 1 > total for w in moduli):
        return False
    if not moduli:
        return True
    modulus = prod(moduli)
    if 27 * modulus**2 > total**3:
        return False
    # At each modulus the exceptional part receives the entire row sum.
    terms = []
    for w in moduli:
        other = modulus // w
        terms.append(total * other * pow(other, -1, w) % modulus)
    # Permute the three slots to put the first exceptional position first.
    for assignment in product(range(3), repeat=len(terms)-1):
        residues = [terms[0], 0, 0]
        for i, term in zip(assignment, terms[1:]):
            residues[i] += term
        minima = [(x-1) % modulus + 1 for x in residues]
        if sum(minima) > total:
            continue
        if bases:
            assert modulus >= total and sum(minima) == total
            parts = tuple(scale * a for a in minima)
            if any(valuation(scale * total, parts, p) for p in bases):
                continue
        return True
    return False


def pair_stage(n, weight, predecessor):
    """0 means no residue pair, 1 needs higher carries, 2 passes the short-row check."""
    best = 0
    weight_factors = list(weight.items())
    if len(weight_factors) == 1:
        return 0
    # A second digit unit in the same place can give a gap of one.
    gap_moduli = []
    for p, power in weight_factors:
        tail = n // power
        w = 1
        if tail % p == 1:
            w = p
            tail //= p
            while tail % p == 0:
                w *= p
                tail //= p
        gap_moduli.append(w)

    @lru_cache(maxsize=None)
    def row(total, moduli, bases=()):
        return permits_row(total, moduli, n // total, bases)

    predecessor_powers = tuple(predecessor.values())
    for tail_flags in product(range(2), repeat=len(weight_factors)-1):
        flags = (0,) + tail_flags
        u = prod(power for (_, power), side in zip(weight_factors, flags) if side == 0)
        v = n // u
        if min(u, v) < 3:
            continue
        gaps_a = tuple(w for w, side in zip(gap_moduli, flags) if w > 1 and side == 0)
        gaps_c = tuple(w for w, side in zip(gap_moduli, flags) if w > 1 and side == 1)
        if not row(v, gaps_a) or not row(u, gaps_c):
            continue
        for places in product(range(2), repeat=len(predecessor_powers)):
            mods_a = gaps_a + tuple(w for w, side in zip(predecessor_powers, places) if side == 0)
            if not row(v, mods_a):
                continue
            mods_c = gaps_c + tuple(w for w, side in zip(predecessor_powers, places) if side == 1)
            if not row(u, mods_c):
                continue
            best = 1
            assert prod(mods_a) >= v or prod(mods_c) >= u
            for total, mods, side in ((v, mods_a, 0), (u, mods_c, 1)):
                if prod(mods) < total:
                    continue
                bases = tuple(p for (p, _), s in zip(weight_factors, flags) if s == side)
                bases += tuple(p for p, s in zip(predecessor, places) if s == side)
                if not row(total, mods, bases):
                    break
            else:
                return 2
    return best


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--start", type=int, default=10**12 + 1234567890)
    parser.add_argument("--count", type=int, default=100000)
    args = parser.parse_args()
    if args.start < 3 or args.count < 1:
        parser.error("need start >= 3 and count >= 1")
    # Examples selected by each of the four sufficient tests in the note.
    for n in (1001234644501, 1001234622970, 1001234650368, 1001234651976):
        m, h = factorizations(n-1, n)
        assert pair_stage(n, h, m) == 0, n
    print("All four examples exclude every normalized row-pair assignment.", flush=True)
    assert permits_row(222788354 // 283, (29, 127, 241), 283, (29, 127, 241, 283))

    factors = factorizations(args.start-1, args.start+args.count-1)
    residue = selected = high = 0
    for i in range(args.count):
        n = args.start + i
        tag, _ = criterion(n, factors[i+1], factors[i])
        stage = pair_stage(n, factors[i+1], factors[i])
        if tag not in ("outside", "inconclusive"):
            assert stage == 0, (n, tag)
            selected += 1
        residue += stage >= 1
        if stage == 2:
            high += 1
            print(f"Higher-carry survivor: {n}", flush=True)
    print(f"{args.count} weights, {selected} sufficient certificates cross-checked, "
          f"{residue} weights with a feasible residue assignment, {high} pass higher carries.")
