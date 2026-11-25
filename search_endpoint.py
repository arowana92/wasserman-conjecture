"""Decide the (1,1,N-2) subcase by splitting prime-power constraints between rows."""

from argparse import ArgumentParser
from itertools import product
from math import comb, gcd, prod, isqrt


def smallest_factors(limit):
    spf = list(range(limit + 1))
    for p in range(2, isqrt(limit) + 1):
        if spf[p] == p:
            for j in range(p*p, limit+1, p):
                if spf[j] == j:
                    spf[j] = p
    return spf


def atoms(n, spf):
    """Return the distinct prime bases together with their full powers in n."""
    result = []
    while n > 1:
        p = spf[n]
        power = 1
        while n % p == 0:
            n //= p
            power *= p
        result.append((p, power))
    return result


def splits(factors):
    for mask in range(1 << len(factors)):
        left = [item for j, item in enumerate(factors) if mask & (1 << j)]
        right = [item for j, item in enumerate(factors) if not mask & (1 << j)]
        yield prod(v for _, v in left), left, right


def avoids(row, p):
    a, b, c = row
    while a or b or c:
        a, x = divmod(a, p)
        b, y = divmod(b, p)
        c, z = divmod(c, p)
        if x+y+z >= p:
            return False
    return True


def find_row(total, scale, scale_atoms, shift, shift_atoms):
    # Look for parts scale*y_i, with sum(y_i) = total/scale.
    target = total // scale
    # At least two y_i must be multiples of each assigned prime power of N-1.
    if any(2*w+1 > target for _, w in shift_atoms) or 27*shift**2 > target**3:
        return None
    bases = [p for p, _ in scale_atoms + shift_atoms]
    # This CRT element is target mod w, and 0 at all other prime powers.
    units = [(target * (shift//w) * pow(shift//w, -1, w)) % shift for _, w in shift_atoms]
    for groups in product(range(3), repeat=len(units)):
        # Number bins by the first atom they receive, allowing unused bins.
        if groups and (groups[0] != 0 or (2 in groups and
                       (1 not in groups or groups.index(2) < groups.index(1)))):
            continue
        residues = [0, 0, 0]
        for j, unit in zip(groups, units):
            residues[j] += unit
        a0, b0, c0 = (value % shift or shift for value in residues)
        left = target - a0 - b0 - c0
        if left < 0:
            continue
        # By construction the residue sum is target mod shift.
        for da in range(0, left+1, shift):
            for db in range(0, left-da+1, shift):
                row = (scale*(a0+da), scale*(b0+db), scale*(c0+left-da-db))
                if all(avoids(row, p) for p in bases):
                    return tuple(sorted(row))
    return None


def candidate(n, spf):
    factors = atoms(n, spf)
    shifted = atoms(n-1, spf)
    if len(shifted) < 3 or len(factors) < 2:
        return None
    heavy_splits = [s for s in splits(shifted) if len(s[1]) >= 3]
    for u, u_factors, v_factors in splits(factors):
        v = n // u
        if min(u, v) < 3:
            continue
        for s, s_factors, t_factors in heavy_splits:
            if s < v:
                continue
            t = (n-1) // s
            first = find_row(n, u, u_factors, s, s_factors)
            if first is None:
                continue
            second = find_row(n, v, v_factors, t, t_factors)
            if second is not None:
                return ((1, 1, n-2), first, second)
    return None


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('first', type=int)
    parser.add_argument('last', type=int)
    args = parser.parse_args()
    if not 3 <= args.first <= args.last:
        parser.error('require 3 <= first <= last')
    spf = smallest_factors(args.last)
    for n in range(args.first, args.last+1):
        rows = candidate(n, spf)
        if rows:
            if any(min(row) < 1 or sum(row) != n for row in rows):
                raise ArithmeticError('Invalid row')
            print('Candidate', n, rows, flush=True)
            values = [comb(n, a)*comb(n-a, b) for a, b, _ in rows]
            if gcd(*values) != 1:
                raise ArithmeticError('Candidate failed the exact coefficient check')
            print('Exact coefficient gcd is one. Independently recheck these rows.')
            return
    print(f'Excluded counterexamples containing (1,1,N-2) for {args.first}..{args.last}.')


if __name__ == '__main__':
    main()
