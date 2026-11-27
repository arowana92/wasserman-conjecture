"""Check the predecessor digit criteria using exact coefficients and polynomial arithmetic."""

from itertools import combinations_with_replacement, permutations, product
from math import comb, gcd, isqrt

from check_endpoint_pair_bound import prime_powers, rows


def digits(n, base):
    result = []
    while n:
        n, digit = divmod(n, base)
        result.append(digit)
    return result


def sparse_profile(total, p, d):
    """Return zero for one unit, the ratio for two units, or None."""
    units = []
    layer = 1
    for digit in digits(total, p) + [0] * len(digits(d, p)):
        if digit:
            if len(units) + digit > 2:
                return None
            units.extend([layer] * digit)
        layer *= p
        if units and d*(total % layer) <= layer:
            if len(units) == 1:
                return 0
            return units[1] // units[0]
    return None


def good_factors(f):
    """Return factors 1+tX with 2 <= t <= f[1] and nonnegative quotient."""
    result = {}
    for t in range(2, f[1]+1):
        if f[-1] % t:
            continue
        quotient = [1]
        for digit in f[1:-1]:
            quotient.append(digit-t*quotient[-1])
        if min(quotient) >= 0 and t*quotient[-1] == f[-1]:
            result[t] = quotient
    return result


def detects(d, p, P, q, R):
    h = (d*R) % P
    ratio = sparse_profile(d*R, p, d)
    other = sparse_profile(d*P, q, d)
    if ratio is not None and (ratio == 0 or other is not None
                              or d <= (1+ratio)*(4*(1+ratio)+1)):
        return True
    if not h or P <= d or d*h > P:
        return False
    units = sum(digits(h, p))
    if gcd(d, h) > 3 and units > 2:
        return False
    f = digits(1+d*P*R, P)
    if d*(h+2)*max(f) <= 4*P and not good_factors(f):
        return True
    m = d*R // P
    if not m or d*d*(h*h+4*m) > 32*P:
        return False
    discriminant = h*h-4*m
    return h <= 3 or units <= 2 or discriminant < 0 or isqrt(discriminant)**2 != discriminant


def exact_triples(bound):
    prime_list = [p for p in range(2, bound+1) if list(prime_powers(p)) == [(p, p)]]
    weights = count = intersections = 0
    for n in range(3, bound+1):
        choices = []
        for p, full in prime_powers(n-1):
            power = p
            while power <= full:
                choices.append((p, power))
                power *= p
        if not any(detects((n-1)//(P*R), p, P, q, R)
                   for (p, P), (q, R) in permutations(choices, 2) if p != q):
            continue
        weights += 1
        minimal = set()
        primes = [r for r in prime_list if r <= n]
        for _, coefficient in rows(n):
            mask = sum(1 << i for i, r in enumerate(primes) if coefficient % r == 0)
            count += 1
            if not any(old & mask == old for old in minimal):
                minimal = {old for old in minimal if old & mask != mask}
                minimal.add(mask)
        for a, b, c in combinations_with_replacement(minimal, 3):
            assert a & b & c, n
            intersections += 1
    print(f"Exact three-way gcd check: {weights} weights, {count} coefficients, "
          f"{intersections} minimal-support triples through {bound}.")


def grouped_cuts():
    cases = allowed = 0
    for p in (2, 3, 5, 7):
        for layer in (p, p*p, p**3, p**4):
            for d in range(2, 27):
                for R in range(2, 64):
                    h = d*R % layer
                    if R % p == 0 or not h or d*h > layer:
                        continue
                    for a in range(1, d):
                        cases += 1
                        if comb(d*R, a*R) % p == 0:
                            continue
                        allowed += 1
                        v = a*R % layer
                        assert 0 < v < h and d*v == a*h
                        assert a % (d//gcd(d, h)) == 0
                        if sum(digits(h, p)) <= 2:
                            assert sum(digits(v, p)) == 1
                            assert next(iter(prime_powers(v)), (p, 1))[0] == p
    assert allowed
    print(f"Grouped projection: {allowed} carry-free cuts among {cases} exact binomials.")


def polynomial_controls():
    checked = surviving = 0
    # Independent of division: generate all eligible nonnegative products,
    # then compare complete coefficient lists against the division test.
    for degree in range(2, 6):
        generated = {}
        for t in range(2, 6):
            for tail in product(range(6), repeat=degree-1):
                g = (1, *tail)
                f = tuple(a+t*b for a, b in zip((*g, 0), (0, *g)))
                if not f[-1] or max(f) > 5:
                    continue
                generated.setdefault(f, {})[t] = list(g)
        for middle in product(range(6), repeat=degree-1):
            for lead in range(1, 6):
                f = (1, *middle, lead)
                candidates = good_factors(f)
                assert candidates == generated.get(f, {}), f
                checked += 1
                surviving += bool(candidates)
    assert surviving
    print(f"Polynomial controls: {checked} digit polynomials, {surviving} admissible factorizations by weight.")


def positive_split():
    # This pair takes the two predecessor bases separately. Coverage of
    # the other primes is deliberately not required by the control.
    P, R, d = 347, 17, 41
    U, V = P+1, 2*P+1
    n = U*V
    A = (P*U, P*U, U)
    C = (5*R*V, 15*R*V, n-20*R*V)
    assert n == 1+d*P*R and gcd(U, V) == 1
    for base, row in ((P, A), (R, C)):
        assert min(row) > 0 and sum(row) == n
        assert all(sum(a % layer for a in row) == n % layer for layer in prime_layers(base, n))
    f = digits(n, P)
    assert d*(f[1]+2)*max(f) <= 4*P and good_factors(f) == {2: [1, 1]}
    assert A[0] + A[2] == U*U and A[2] == U
    print(f"Positive split control at {n} with coprime scales {U}, {V}.")


def prime_layers(p, n):
    w = p
    while w <= n:
        yield w
        w *= p


def examples():
    assert pow(41, 1024, 8192) == 1 and 385*41**2 % 8192 == 17
    assert sparse_profile(385*41**2, 2, 385) == 16
    P, R, d = 1024, 20610923369137, 55
    assert list(prime_powers(R)) == [(R, R)]
    n = 1+d*P*R
    f = digits(n, P)
    assert f == [1, 7, 1, 2, 4, 7, 1] and detects(d, 2, P, R, R)
    assert sparse_profile(d*R, 2, d) is None and sparse_profile(d*P, R, d) is None
    assert d*(f[1]+2)*max(f) <= 4*P and gcd(d, f[1]) == 1
    print(f"Five-base exponent family and degree-six example at {n} verified.")


if __name__ == "__main__":
    grouped_cuts()
    polynomial_controls()
    positive_split()
    examples()
    exact_triples(250)
