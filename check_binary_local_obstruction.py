"""Construct simultaneous low-digit witnesses on the cofactor 12031."""
from fractions import Fraction
from math import gcd, lcm

from prime_certificates import primes_through


D = 12031
SLOPES = (681, 9222, 2128)
SEED = 106675
PRIME_SEED = 265327
STEP = 22211280
TARGETS = {3: 17140160990152208, 53: 33, 227: 8}
SKIP = (5, 7, 13, 113, 1093, 3511)


def row(q):
    x, y, z = SLOPES
    return x*q, y*q, z*q+1


def valuation(parts, p):
    total = sum(parts)
    w = p
    v = 0
    while w <= total:
        v += total//w - sum(k//w for k in parts)
        w *= p
    return v


def check_seeds():
    x, y, z = SLOPES
    assert D == 53*227 == x+y+z and (x&y | y&z | z&x) == 0
    assert STEP == 16*27*5*7*13*113
    assert SEED % 5876 == 907 and STEP % 5876 == 0
    for r in SKIP:
        assert STEP % (r-1) == 0 and (D*pow(2, SEED, r)+1) % r
    for p, q in TARGETS.items():
        m = p
        while STEP % m == 0:
            m *= p
        assert pow(2, SEED, m) == q % m
        assert pow(2, STEP, m) == 1 and pow(2, STEP, m*p) != 1
        assert valuation(row(q), p) == 0
    assert (D*TARGETS[3]+1) % 9 == 6
    assert [r for r in primes_through(10**6) if r > 2 and pow(2,r-1,r*r) == 1] == [1093,3511]


def lift(a, period, p, target, depth):
    """Impose 2^a = target modulo p^depth while preserving the old progression."""
    modulus = p
    while modulus < p**depth and pow(2, period, modulus*p) == 1:
        modulus *= p
    assert pow(2, a, modulus) == target % modulus and pow(2, period, modulus) == 1
    while modulus < p**depth:
        next_modulus = p*modulus
        value = pow(2, a, next_modulus)
        derivative = (pow(2, period, next_modulus)-1)//modulus
        digit = ((target-value)//modulus)*pow(value*derivative, -1, p) % p
        a += digit*period
        period *= p
        modulus = next_modulus
    return a, period


def witness(primes, depth, first_layers=None, extra=1):
    a, period = (PRIME_SEED if first_layers else SEED), STEP
    for p, q in TARGETS.items():
        precision = depth+1
        if first_layers:
            first = first_layers[p]
            assert first >= max(45, depth+3) and extra >= 1
            low = first-(2 if p == 3 else 1)
            q -= p**low
            precision = first+8+extra
        a, period = lift(a, period, p, q, precision)
        if first_layers:
            assert gcd(a, period) == 1
    for r in sorted(set(primes) - {2, *TARGETS, *SKIP}):
        if pow(2, r-1, r*r) == 1 or period % r == 0:
            raise ValueError(("prime outside the induction hypotheses", r))
        new_period = lcm(period, r-1)
        if first_layers:
            # Choose a unit above the previous unit class, even at new factors.
            free = new_period
            while gcd(free, period) > 1:
                free //= gcd(free, period)
            a += period*((1-a)*pow(period, -1, free) % free)
        period = new_period
        if (D*pow(2,a,r)+1) % r == 0:
            modulus = r**(depth+1)
            for k in ((1,2) if first_layers else (1,)):
                target = (3*k*r-1)*pow(D, -1, modulus) % modulus
                candidate, stride = lift(a, period, r, target, depth+1)
                if not first_layers or candidate % r:
                    a, period = candidate, stride
                    break
            else:
                raise ArithmeticError(("no unit lift", r))
        if first_layers:
            assert gcd(a, period) == 1
    return a, period


def check_witness(e, primes, depth):
    # Work directly with the actual row residues, including C via division by 3.
    active = []
    for p in primes:
        modulus = p**(depth+1)
        q = pow(2, e, 3*modulus)
        n = D*q+1
        assert n % 3 == 0
        if p > 2 and n % p and (n-1) % p:
            continue  # The full endpoint coefficient is a unit.
        if p != 2 and p not in TARGETS:
            assert n % (p*p) in {3*p % (p*p), 6*p % (p*p)}
            active.append(p)
        parts = (1, 1, n-2), row(q), (n//3,)*3
        assert any(all(sum(k % (p**i) for k in triple) < p**i
                           for i in range(1, depth+2)) for triple in parts)
    return active


def check_prime_seed():
    assert gcd(PRIME_SEED,STEP) == 1 and PRIME_SEED % 5876 == 907
    assert [(D*pow(2,PRIME_SEED,r)+1) % r for r in SKIP] == [4, 4, 2, 84, 1085, 3254]
    for (p, q), residue in zip(TARGETS.items(), (1,31,3)):
        a, h = lift(PRIME_SEED, STEP, p, q, 6)
        assert h % p == 0 and a % p == residue


def check_middle(e, first_layers, extra):
    for p, first in first_layers.items():
        modulus = p**(first+8+extra)
        parts = row(pow(2,e,modulus))
        layers = []
        w = p
        while w <= modulus:
            layers.append(sum(x % w for x in parts)//w)
            w *= p
        assert next(i for i, cost in enumerate(layers, 1) if cost) == first
        assert all(layers[i-1] == 2 for i in range(first+8,first+8+extra))


def leading_ternary_costs(e, count):
    """Enclose 2^e by integer mantissas times a common power of three."""
    precision = e.bit_length()+24
    block = 3**precision
    lower = upper = 1
    shift = 0
    for bit in bin(e)[2:]:
        factor = 2 if bit == '1' else 1
        lower, upper = factor*lower**2, factor*upper**2
        shift *= 2
        while upper >= block**2:
            lower //= block
            upper = (upper+block-1)//block
            shift += precision

    # At the top layer 3^(shift+j), determine j by integer comparison.
    lo, hi = 0, 2*precision+20
    assert D*lower >= 1 and D*upper+1 < 3**hi
    while hi-lo > 1:
        mid = (lo+hi)//2
        if 3**mid <= D*lower:
            lo = mid
        else:
            hi = mid
    w = 3**lo
    assert D*upper+1 < 3*w and lo >= count-1
    costs = []
    for _ in range(count):
        floors = []
        for slope, offset in zip((*SLOPES,D), (0,0,1,1)):
            # The normalized offset is offset / 3^shift and lies in [0,1].
            left = (slope*lower+(offset if shift == 0 else 0))//w
            right = (slope*upper+offset)//w
            assert left == right
            floors.append(left)
        costs.append(floors[-1]-sum(floors[:-1]))
        w //= 3
    return costs


def check_explicit_prime():
    for exponent in (200,1000,2000):
        parts = row(2**exponent)
        n = sum(parts)
        w = 1
        while w*3 <= n:
            w *= 3
        costs = []
        for _ in range(4):
            costs.append(n//w-sum(x//w for x in parts))
            w //= 3
        assert leading_ternary_costs(exponent,4) == costs

    primes = primes_through(19)
    layers = dict.fromkeys(TARGETS,45)
    a, h = witness(primes, 2, layers, 1)
    factor = 11**246
    b = ((a-1)*pow(factor,-1,h)) % h
    e = 1+factor*(b+4430*h)
    assert e % h == a % h and e.bit_length() == 1694
    assert (e-1) % factor == 0 and factor**2 > e
    assert pow(2,e-1,e) == 1 and gcd(pow(2,(e-1)//11,e)-1,e) == 1
    check_witness(e, primes, 2)
    check_middle(e, layers, 1)
    assert leading_ternary_costs(e,4) == [0,0,0,0]
    return e


def leading_windows():
    # The floor difference is periodic modulo one. Check the whole rational
    # forward orbit, which certifies windows of every finite leading length.
    witnesses = ((3, Fraction(423213851608691,2701703435345984178), 39),
                 (53, Fraction(21,146068), 3), (227, Fraction(3,226), 1))
    for p, t, orbit_size in witnesses:
        assert 1 < D*t < p
        seen = set()
        y = t
        while y not in seen:
            assert (D*y)//1 == sum((s*y)//1 for s in SLOPES)
            seen.add(y)
            y = (p*y) % 1
        assert len(seen) == orbit_size


def first_full_predecessor():
    # At this exponent the odd row covers every prescribed lowest layer.
    q, n = 2**119, 1903*2**119+1
    parts = (264*q, 519*q, 1120*q+1)
    assert sum(parts) == n and n % 9 == 3
    assert valuation(parts, 2) == 0 and valuation(parts, 3) == 81
    assert all(sorted(v % w for v in parts) == [0, 0, 1] for w in (q, 11, 173))
    assert all(v % 3 == 0 for v in parts)
    assert sum(v % 9 for v in parts) >= 9
    assert 1903 == 11*173 and all(pow(2,2580,w) == 1 for w in (9,11,173))


if __name__ == "__main__":
    check_seeds()
    first_full_predecessor()
    leading_windows()
    for bound, depth in [(200, 32), (4096, 12)]:
        primes = list(primes_through(bound))
        a, h = witness(primes, depth)
        for e in (a, a+h):
            active = check_witness(e, primes, depth)
        print(f"Primes through {bound}, depth {depth} beyond the first layer, "
              f"period of {h.bit_length()} bits, additional weight primes {active}.")

    check_prime_seed()
    for bound, depth, first, extra in [(19, 2, (45,45,45), 1),
                                       (4096, 12, (60,53,71), 8)]:
        primes = primes_through(bound)
        layers = dict(zip(TARGETS,first))
        a, h = witness(primes, depth, layers, extra)
        for e in (a,a+h):
            check_witness(e, primes, depth)
            check_middle(e, layers, extra)
        print(f"Primitive progression through prime {bound}, first layers {first}, "
              f"{extra} prescribed double layers per base, period of {h.bit_length()} bits.")

    e = check_explicit_prime()
    print(f"Certified prime exponent of {e.bit_length()} bits, first layers (45,45,45), "
          "four carry-free highest ternary layers.")
