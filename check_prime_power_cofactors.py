"""Eliminate endpoint counterexamples at N=dQ+1 using finite carry automata."""

from argparse import ArgumentParser
from collections import deque
from fractions import Fraction
from math import comb, gcd, prod

from search_endpoint import atoms, avoids, smallest_factors


def shapes(d, factors):
    """Generate all normalized row pairs meeting the lowest-layer conditions."""
    for u in range(3, (d - 1) // 2 + 1):
        if gcd(u, d) != 1:
            continue
        for a in range(u, (d - 1) // 2 + 1, u):
            for b in range(a, d - a, u):
                # A=(aQ,bQ,cQ+1). If just one unshifted part is zero
                # modulo w, the other must be a unit to have residue one.
                required = [w for p, w in factors if not (
                    (a % w == 0 and (b % w == 0 or b % p != 0)) or
                    (b % w == 0 and a % p != 0))]
                if any(2*w + 1 > u for w in required) or 27*prod(required)**2 > u**3:
                    continue
                for f in range(1, u // 3 + 1):
                    for g in range(f, (u - f) // 2 + 1):
                        h = u - f - g
                        if all(sum(k % w == 0 for k in (f, g, h)) == 2 for w in required):
                            yield u, a, b, d-a-b, f, g, h


def carry_language(p, slopes, initial, first_digits):
    """Does some allowed nonnegative input yield carry-free affine parts?"""
    queue = deque()
    seen = set()

    def advance(state, digits):
        for digit in digits:
            values = tuple(s + k*digit for s, k in zip(state, slopes))
            if sum(v % p for v in values) >= p:
                continue
            next_state = tuple(v // p for v in values)
            if next_state not in seen:
                seen.add(next_state)
                queue.append(next_state)

    advance(initial, first_digits)
    while queue:
        state = queue.popleft()
        if avoids(state, p):
            return True
        advance(state, range(p))
    return False


def locally_possible(d, shape, p, binary_q, residue, stride):
    u, a, b, c, f, g, h = shape
    if p == 2 and binary_q:
        # Q=2^e shifts the three multipliers, and +1 occupies a fresh bit.
        return u % 2 == 1 and not ((a & b) | (a & c) | (b & c))
    v0 = (d*residue + 1) // u
    v_step = d*(stride // u)
    digits = [x for x in range(p)
              if ((p != 2 and d % p and not binary_q) or (residue + stride*x) % p)
              and (u % p or (v0 + v_step*x) % p)]
    if u % p == 0 or d % p == 0:
        if carry_language(p, (a*stride, b*stride, c*stride),
                          (a*residue, b*residue, c*residue+1), digits):
            return True
    if u % p:
        return carry_language(p, (f*v_step, g*v_step, h*v_step),
                              (f*v0, g*v0, h*v0), digits)
    return False


def survives(d, factors, shape, spf, binary_q):
    u, a, b, _, f, g, h = shape
    q0 = -pow(d, -1, u) % u
    bases = sorted({2} | {p for p, w in factors + atoms(u, spf)})
    if not all(locally_possible(d, shape, p, binary_q, q0, u) for p in bases):
        return False
    if not binary_q:
        # When the other row cannot handle w, a single unshifted
        # nonzero residue of A has to equal 1. Impose it in base 2 too.
        for p, w in factors:
            if p == 2 or sum(k % w == 0 for k in (f, g, h)) >= 2:
                continue
            if a % w == b % w == 0:
                continue
            q1 = pow(a+b, -1, w)
            residue = q0 + u*((q1 - q0)*pow(u, -1, w) % w)
            if not locally_possible(d, shape, 2, False, residue, u*w):
                return False
    return True


def check_automaton():
    # All integers t yield a fixed row when the slopes vanish.
    for n in range(3, 31):
        for a in range(1, n//3 + 1):
            for b in range(a, (n-a)//2 + 1):
                row = a, b, n-a-b
                coefficient = comb(n, a)*comb(n-a, b)
                for p in (2, 3, 5):
                    assert carry_language(p, (0, 0, 0), row, range(p)) == (coefficient % p != 0)
    # A carry-free row (x,x,1) in base 2 requires x=0 for any input.
    assert not carry_language(2, (1, 1, 0), (1, 1, 1), range(2))
    # t=0 is rejected for x=t+3; the even input t=6 gives 9,18,36.
    assert not avoids((3, 6, 12), 2)
    assert carry_language(2, (1, 2, 4), (3, 6, 12), [0])
    for delta in (0, 1):
        exponents = {}
        for h in range(12):
            if h == 1:
                continue
            for j in range(h+1):
                slope = Fraction(j-delta, h-1)
                exponents[slope] = max(12-h, exponents.get(slope, 0))
        assert len(exponents) == 43 and sum(exponents.values()) == 202


def main():
    arg = ArgumentParser(description=__doc__)
    arg.add_argument("limit", nargs="?", type=int, default=1000)
    limit = arg.parse_args().limit
    check_automaton()
    spf = smallest_factors(limit)
    total = 0
    survivors = {False: set(), True: set()}
    for d in range(1, limit+1):
        factors = atoms(d, spf)
        for shape in shapes(d, factors):
            total += 1
            for binary in ([False, True] if d % 2 else [False]):
                if survives(d, factors, shape, spf, binary):
                    survivors[binary].add(d)
                    if d % 2:
                        print(f"Surviving odd cofactor shape: d={d}, {shape}, Q even={binary}")
    print(f"Checked {total} normalized row pairs for 1 <= d <= {limit}.")
    print("Possible cofactors for odd Q:", sorted(survivors[False]))
    print("Possible cofactors for Q a power of two:", sorted(survivors[True]))


if __name__ == "__main__":
    main()
