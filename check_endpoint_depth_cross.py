"""Independent cut enumeration and factorial-floor checks for truncated Kummer tests."""

from itertools import combinations_with_replacement
from math import comb, prod

from search_endpoint import atoms, smallest_factors, splits


RECORDS = (
    (771, (210, 231, 330), (257, 257, 257), 2310),
    (97876, (14500, 27000, 56376), (24469, 24469, 48938), 391500),
    (2556249, (576264, 958593, 1021392), (852083,) * 3, 24),
    (25008076, (5594400, 7409800, 12003876), (6252019, 6252019, 12504038), 2700),
)


def floor_cost(n, row, q):
    return n // q - sum(x // q for x in row)


def trial(n):
    p = 2
    while p * p <= n:
        power = 1
        while n % p == 0:
            n //= p
            power *= p
        if power > 1:
            yield p, power
        p = 3 if p == 2 else p + 2
    if n > 1:
        yield n, n


def check_record(n, a, b, gcd_expected, depth):
    assert sum(a) == sum(b) == n and min(*a, *b) > 0
    witnesses = []
    endpoint_gcd = 1
    for p, power in (*trial(n), *trial(n - 1)):
        first = []
        valuations = []
        for row in (a, b):
            q, cost, initial = p, 0, None
            while q <= n:
                v = floor_cost(n, row, q)
                if v and initial is None:
                    initial = q
                cost += v
                q *= p
            first.append(initial)
            valuations.append(cost)
        if None not in first:
            q, extra = power, 0
            while q < max(first):
                q *= p
                extra += 1
            witnesses.append((p, extra))
            endpoint_gcd *= min(power, p**min(valuations))
    assert min(extra for _, extra in witnesses) == depth
    assert endpoint_gcd == gcd_expected
    return witnesses


def admits_row(n, u, weight, pred, depth):
    shift = prod(q for _, q in pred)
    cap = n // u
    if cap < 3 or any(2*q >= cap for _, q in pred) or 27 * shift**2 > cap**3:
        return False
    # Enumerate individual actual cuts mod u*shift. Two possible residues per new power.
    residues = [0]
    period = u
    for _, q in pred:
        inv = pow(period, -1, q)
        residues = [z + period * (((v-z)*inv) % q) for z in residues for v in (0, 1)]
        period *= q
    cuts = sorted(a for z in residues for a in range(z or period, n, period))
    allowed = set(cuts)
    levels = []
    for p, power in weight + pred:
        q = p
        while q <= min(n, power * p**depth):
            levels.append(q)
            q *= p
    for a, b in combinations_with_replacement(cuts, 2):
        c = n - a - b
        if c >= b and c in allowed and all(floor_cost(n, (a,b,c), q) == 0 for q in levels):
            return True
    return False


def check_first_weights(limit=100000):
    spf = smallest_factors(limit)
    depth = 0
    first = []
    for n in range(3, limit + 1):
        nf, mf = atoms(n, spf), atoms(n-1, spf)
        if len(nf) < 2 or len(mf) < 3:
            continue
        for u, left, right in splits(nf):
            if u > n // u or min(u, n // u) < 3:
                continue
            for s, picked, rest in splits(mf):
                sides = [(u, left, picked), (n//u, right, rest)]
                # Try the larger actual-part modulus first.
                if u*s < (n//u)*((n-1)//s):
                    sides.reverse()
                def both(k):
                    return all(admits_row(n, *side, k) for side in sides)
                while both(depth):
                    first.append(n)
                    depth += 1
    assert first == [n for n, *_ in RECORDS if n <= limit]
    return first


def analytic_margins():
    # Taylor lower bound for exp(2.54), and e < 68/25.
    from fractions import Fraction
    from math import factorial
    assert sum(Fraction(127, 50)**i / factorial(i) for i in range(12)) > Fraction(25, 2)
    assert Fraction(68, 25) * Fraction(201, 1000) < Fraction(11, 20)
    assert 7**89 * 11**(127*89) * 68**6500 < 20**(127*89) * 25**6500
    assert Fraction(12) + (Fraction(115, 3) - 65)/2 == -Fraction(4, 3)
    assert 18 + 15*Fraction(18, 100) < 21


if __name__ == '__main__':
    analytic_margins()
    for depth, record in enumerate(RECORDS, 1):
        print(record[0], check_record(*record, depth))
    a, b = RECORDS[0][1:3]
    from math import gcd
    assert gcd(771*770, *(comb(771,r[0])*comb(771-r[0],r[1]) for r in (a,b))) == 2310
    print('Independent first weights through 100000', check_first_weights())
