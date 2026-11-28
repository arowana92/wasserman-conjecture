"""Enumerate variable cofactors at the two sparse predecessor layers."""

from bisect import bisect_left, bisect_right
from functools import lru_cache
from itertools import product
from math import gcd

from check_residual_power_classes import paired_binary_ternary_cofactors
from check_three_unit_layers import grouped_table


@lru_cache(maxsize=4)
def powers_through(base, limit):
    powers = [1]
    while powers[-1] < limit:
        powers.append(powers[-1]*base)
    return powers


def exponent_lookup(p, q, max_exponent):
    """For E*p^a = 1 (mod q^h), return candidate a up to the requested bound."""
    tables = []
    modulus = 1
    while True:
        modulus *= q
        inverse = pow(p, -1, modulus)
        residue = 1
        table = {}
        for a in range(1, max_exponent+1):
            residue = residue*inverse % modulus
            table[residue] = a
            if residue == 1:
                break
        tables.append((modulus, table, a))
        if a == max_exponent:
            return tables


def possible_cofactors(inverse, p, margin, upper, powers):
    # Use intervals [B^(j-1), B^j) for B the largest p-power <= margin.
    # At most one cofactor in each interval has the required inverse residue.
    if p == 2:
        # A prefix E with top bit h needs a zero run just above that bit.
        end = (margin*upper-1).bit_length()
        low = inverse & ((1 << end)-1)
        heads = low
        for gap in range(1, (margin-1).bit_length()):
            heads &= ~(low >> gap)
        while heads:
            leading_bit = heads & -heads
            heads -= leading_bit
            E = low & (2*leading_bit-1)
            if E <= upper:
                modulus = 1 << (margin*E-1).bit_length()
                if low & (modulus-1) == E:
                    yield E
        return

    block = powers[bisect_right(powers, margin)-1]
    assert block > 1
    ceiling = block
    while ceiling//block <= upper:
        E = inverse % ceiling
        if ceiling//block <= E <= upper:
            layer = powers[bisect_left(powers, margin*E)]
            if inverse % layer == E:
                yield E
        ceiling *= block


def paired_weights(p, q, s, limit, cofactor_classes=None):
    """All E,a,b with two common-remainder layers and 1+s*E*p^a*q^b <= limit."""
    margin = s*s
    powers = powers_through(p, margin*limit)
    other = powers_through(q, margin*limit)
    max_a = bisect_right(powers, (limit-1)//(s*q))-1
    if max_a < 1:
        return set()
    tables = exponent_lookup(p, q, max_a)
    # Each needed modulus divides this one. Divide the inverse by q exactly
    # after adding the appropriate multiple of the modulus.
    modulus = powers[bisect_left(powers, margin*((limit-1)//(s*p*q)))]
    correction = -pow(modulus, -1, q)
    inverse = 1
    result = set()
    for b in range(1, len(other)):
        Q = other[b]
        upper = (limit-1)//(s*p*Q)
        if upper == 0:
            break
        inverse = (inverse+(inverse*correction % q)*modulus)//q
        for E in possible_cofactors(inverse, p, margin, upper, powers):
            if cofactor_classes is not None and E % 1320 not in cofactor_classes:
                continue
            if E % q == 0:
                continue
            h = bisect_left(other, margin*E)
            layer = other[h]
            small, table, period = tables[min(h,len(tables))-1]
            first = table.get(E % small)
            if first is None:
                continue
            last = bisect_right(powers, (limit-1)//(s*E*Q))-1
            for a in range(first, last+1, period):
                if E*powers[a] % layer == 1:
                    result.add((E,a,b))
    return result


def first_carry(x, y, prime):
    place = prime
    while x or y:
        x, d = divmod(x, prime)
        y, e = divmod(y, prime)
        if d+e >= prime:
            return place
        place *= prime
    return None


def check_table_case(p, q, w, limit):
    s = w+1
    rows = paired_weights(p, q, s, limit)
    assert rows == {(E,a,b) for E,b,a in paired_weights(q, p, s, limit)}
    with_unitary = 0
    remaining = []
    divisors = [u for u in range(3,w+1) if w % u == 0]
    for E,a,b in rows:
        W = E*p**a*q**b
        N = s*W+1
        eligible = [u for u in divisors if (W+1) % u == 0 and N % u == 0 and gcd(u,N//u) == 1]
        if not eligible:
            continue
        with_unitary += 1
        carry_at_p = first_carry(w*E*q**b, E*q**b, p) is not None
        carry_at_q = first_carry(w*E*p**a, E*p**a, q) is not None
        assert carry_at_p == (binomial_valuation(w*W, W+1, p) > 0)
        assert carry_at_q == (binomial_valuation(w*W, W+1, q) > 0)
        if carry_at_p or carry_at_q:
            continue
        remaining.append((E,a,b))
    return len(rows), with_unitary, remaining


def direct_small_check(p, q, s, limit):
    rows = set()
    for a,b in product(range(1,limit.bit_length()), repeat=2):
        P, Q = p**a, q**b
        for E in range(1,(limit-1)//(s*P*Q)+1):
            if gcd(E,p*q) != 1:
                continue
            x, y = E*Q-1, E*P-1
            L, M = 1, 1
            while x % p == 0:
                x //= p
                L *= p
            while y % q == 0:
                y //= q
                M *= q
            if min(L,M) >= s*s*E:
                rows.add((E,a,b))
    assert rows == paired_weights(p,q,s,limit)
    return len(rows)


def binomial_valuation(x, y, prime):
    n = x+y
    total = 0
    while n:
        n //= prime
        x //= prime
        y //= prime
        total += n-x-y
    return total


def extended_eleven_range(limit):
    classes = paired_binary_ternary_cofactors()
    rows = paired_weights(2, 3, 11, limit, classes)
    survivors = []
    for E,a,b in sorted(rows):
        W = (E*3**b) << a
        if W % 25 != 4 or W % 11 != 10:
            continue
        values = [binomial_valuation(10*W, W+1, r) for r in (2,3,5,11)]
        assert min(values) > 0
        survivors.append((E,a,b,values))
    print(len(rows), 'tuples in the 40 cofactor classes;', survivors)
    if limit == 10**10000:
        assert len(rows) == 71
        assert survivors == [
            (43,1976,5473,[4423,619,2252,1545]),
            (43,1976,15713,[12540,619,5851,3892]),
        ]


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('digits', type=int, nargs='?', default=1000)
    parser.add_argument('--binary-ternary', action='store_true',
                        help='use the 40 possible cofactor classes at bases 2 and 3')
    args = parser.parse_args()
    digits = args.digits
    assert digits > 0
    limit = 10**digits
    if args.binary_ternary:
        extended_eleven_range(limit)
        raise SystemExit()
    assert direct_small_check(2,3,2,500000) > 0
    for p,q,s in [(2,3,5),(3,5,4),(5,7,3)]:
        direct_small_check(p,q,s,500000)
    total = 0
    survivors_total = []
    for p,q,w in grouped_table():
        count, with_unitary, survivors = check_table_case(p,q,w,limit)
        total += count
        survivors_total.extend((p,q,w+1,*row) for row in survivors)
        print(p,q,w+1,count,with_unitary,survivors,flush=True)
    print('Total:',total)
    assert not survivors_total
    if digits == 1000:
        assert total == 2245
