"""Check the local exclusions by forward search and exhibit a delayed common carry."""

from fractions import Fraction
from math import gcd, lcm

from flint import fmpz

from check_cofactor_power_indices import (
    EXPECTED, EXTRA, common_indices, local_obstructions, minimal, records,
)
from check_prime_power_cofactors import carry_language
from prime_certificates import primes_through
from search_endpoint import avoids


def check_local(d, shape, prog, r):
    u, x, y, z, *_ = shape
    R, M = prog
    excluded, witness, failure = local_obstructions(d, shape, prog, r)
    assert avoids((witness*x, witness*y, witness*z+1), r)
    base_part = 1
    while M % (r*base_part) == 0:
        base_part *= r
    assert (witness-R) % base_part == 0
    if u == r:
        assert (d*witness+1) % (r*r) != 0
    assert u == r or u % r
    # The failed precision has a unique possible power residue. Test all its roots.
    expected = set()
    if failure is not None:
        root, w = failure
        assert {pow(v, w//r, w) for v in range(R % r, w, r)} == {root}
        precision = max(w, base_part, r*r if u == r else 1)
        for v in range(root, precision, w):
            if (v-R) % base_part or (r == u and (d*v+1) % (r*r) == 0):
                continue
            row = (x*v, y*v, z*v+1)
            layer = r
            while layer <= precision and sum(k % layer for k in row) < layer:
                layer *= r
            if layer > precision:
                assert not carry_language(r, (x, y, z), tuple(k//precision for k in row), range(r))
        expected.add(w//r)
        assert pow(witness, r-1, w//r) == 1
    else:
        assert witness == 1
    for k in range(2, r):
        if all(pow(t, k, r) != R % r for t in range(1, r)):
            expected.add(k)
    assert minimal(expected) == minimal(excluded)


def prime_root(target, e, r, k):
    """Lift every candidate root one digit at a time, including singular roots."""
    roots = [0]
    w = 1
    for _ in range(k):
        roots = [v+i*w for v in roots for i in range(r)
                 if pow(v+i*w, e, r*w) == target % (r*w)]
        w *= r
    assert roots
    return min(roots), w


def first_carry(parts, r):
    w = r
    j = 1
    while w <= sum(parts):
        if sum(parts)//w > sum(x//w for x in parts):
            return j
        w *= r
        j += 1
    return None


def delayed_example():
    # e=52 makes 53 a compulsory weight prime while leaving its higher digits free.
    e, depth = 52, 12
    R, M = 0, 1
    for r in primes_through(67):
        if r in (2, 7, 17, 53, 67):
            k = depth + (r == 17)
            w = r**k
            if r == 2:
                target = 227*pow(1139, -1, w)   # gives V = 65 locally
            elif r == 53:
                target = (7*53-1)*pow(2278, -1, w)
            else:
                target = {7: 4685, 17: 20557 + 17**depth, 67: 238}[r]
            v, w = prime_root(target, e, r, k)
        else:
            v, w = next(i for i in range(1, r) if (2278*pow(i, e, r)+1) % r), r
        R += M*((v-R)*pow(M, -1, w) % w)
        M *= w
    q = R
    while not fmpz(q).is_prime():
        q += M
    P = q**e
    N = 2278*P+1
    assert P % 15946 == 5665 and N % 49 != 0
    V = N//7
    A, C = (833*P, 1407*P, 38*P+1), (V, 2*V, 4*V)
    assert sum(A) == sum(C) == N
    endpoint, common = [], []
    for r in primes_through(67):
        if N*(N-1) % r:
            continue
        endpoint.append(r)
        first, second = first_carry(A, r), first_carry(C, r)
        assert first is None or second is None or max(first, second) > depth
        if first is not None and second is not None:
            common.append(r)
    assert endpoint == common == [2, 7, 17, 53, 67]
    assert first_carry(A, 17) == 13
    print(f"Prime base {q} at exponent 52. Small common primes {common}; "
          f"every one needs a carry past layer {depth} in at least one row.")
    print(f"A first carries at layer 13 at 17. Weight bit length {N.bit_length()}.")


def verify():
    small, extra = [], []
    count = 0
    for d, shape, prog, first in records():
        barred = {2} if d <= 1000 else {2756}
        for r in first:
            if r == 2:
                continue
            check_local(d, shape, prog, r)
            barred.update(local_obstructions(d, shape, prog, r)[0])
            count += 1
        (small if d <= 1000 else extra).append(minimal(barred))
    assert common_indices(small) == EXPECTED and common_indices(extra) == EXTRA
    # Count exact divisibility patterns among divisors of the common period.
    for indices, target in [(EXPECTED, Fraction(299026867069663408869, 561053121385651819925)),
                            (EXTRA, Fraction(287087939, 653714243))]:
        period = lcm(*indices)
        divisors = [1]
        for p, j in fmpz(period).factor():
            divisors = [x*int(p)**k for x in divisors for k in range(j+1)]
        total = sum(int(fmpz(period//t).euler_phi()) for t in divisors
                    if any(t % n == 0 for n in indices))
        assert Fraction(total, period) == target
    assert all(pow(x, 2756, 53**2) == 1 for x in range(1, 53**2) if x % 53)
    assert (2279//53)*pow(7, -1, 53) % 53 == 44
    assert all(not avoids(s, 53) for s in ((833, 1407, 38), (938, 1309, 31)))
    assert avoids((224*147, 455*147, 101*147+1), 2) and 147 % 4 == 3
    print(f"Checked {count} local obstructions and witnesses, both densities, and the prime 53.")


if __name__ == "__main__":
    verify()
    delayed_example()
