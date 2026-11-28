"""Bounded prime bases, unbounded binary prefixes and exponents up to 10**100."""

from math import gcd, isqrt

from check_even_cofactor_partners import CASES, FixedPrefix, order
from check_even_two_power_cofactors import merge_progressions, prefix_graph


def logarithm(q, residue, modulus, group_size):
    period = order(q, modulus, group_size)
    width = isqrt(period) + 1
    first = {}
    x = 1
    for j in range(width):
        first.setdefault(x, j)
        x = x*q % modulus
    back = pow(x, -1, modulus)
    x = residue % modulus
    for i in range(width + 1):
        if x in first:
            value = (i*width + first[x]) % period
            assert pow(q, value, modulus) == residue % modulus
            return value, period
        x = x*back % modulus
    return None


def exponent_class(q, u, P, residue, modulus):
    if gcd(q, modulus) > 1:
        return None
    binary = modulus & -modulus
    odd = modulus // (binary * P)
    combined = (0, 1)
    for base, power in ((P, P), (u, odd), (2, binary)):
        hit = logarithm(q, residue, power, (power//base)*(base-1))
        if hit is None:
            return None
        combined = merge_progressions(*combined, *hit)
        if combined is None:
            return None
    return combined


def no_carry(entries, p):
    power = p
    while power <= sum(entries):
        if sum(x % power for x in entries) >= power:
            return False
        power *= p
    return True


def small_bases():
    primes = [q for q in range(2, 10001) if all(q % j for j in range(2, isqrt(q)+1))]
    cases = {}
    for w, P, u, a, b, r, residue, modulus in CASES:
        D = w*P
        slopes = (a*P, b*w, (w-a)*P-b*w)
        language = FixedPrefix(2, w, r, prefix_graph(2, slopes, w, r))
        counts = [1]*len(language.edges)
        for _ in range(10000):
            counts = [sum(counts[j] for _, j in row) for row in language.edges]
        bits = max(counts).bit_length()
        assert bits == {15424: 1129, 24352: 849, 27904: 424}[D]
        print(D, "10,000-bit block bound below 2 to the power", bits)
        for q in primes:
            if not no_carry(slopes, q):
                continue
            hit = exponent_class(q, u, P, residue, modulus)
            if hit is not None:
                cases[D, q] = hit
    assert cases == {(27904, 379): (2668154853, 4529848320),
                     (27904, 2539): (1366259897, 1509949440),
                     (27904, 2579): (2879938563, 4529848320),
                     (27904, 7043): (2079006999, 4529848320),
                     (27904, 7507): (388436403, 4529848320)}
    t, L = cases.pop((27904, 7043))
    assert pow(7043, L, 81) == 1
    v = ((1 + 27904*pow(7043, t, 81)) * pow(11, -1, 81)) % 81
    assert v == 39
    for f in range(1, 4):
        for g in range(f, (11-f)//2 + 1):
            shape = (f, g, 11-f-g)
            assert any(sum((v*x) % z for x in shape) >= z for z in (3, 9, 27, 81))
    for (_, q), (t, L) in cases.items():
        if q == 2539:
            assert L % 45 == 0 and t % 5 == t % 9 == 2
            assert 2 not in {i*i % 5 for i in range(5)}
            assert 2 not in {i**3 % 9 for i in range(9)}
        else:
            assert L % 9 == 0 and gcd(t, 9) == 3
    return cases


def lift(q, e, period):
    w, P, u, a, b, r, residue, modulus = CASES[-1]
    slopes = (a*P, b*w, (w-a)*P-b*w)
    edges = FixedPrefix(2, w, r, prefix_graph(2, slopes, w, r)).edges
    transitions = [dict(row) for row in edges]
    binary = modulus & -modulus
    state = 0
    for j in range(binary.bit_length()-1):
        state = transitions[state][(residue >> j) & 1]
    assert 0 < e < period and pow(q, e, binary) == residue % binary
    assert pow(q, period, 2*binary) == binary+1
    bound = 10**100
    strides = []
    stride = period
    while stride <= bound:
        strides.append(stride)
        stride *= 2
    high = binary << len(strides)
    jumps = [pow(q, stride, high) for stride in strides]
    stack = [(0, state, e, pow(q, e, high))]
    terminals = 0
    while stack:
        h, state, n, z = stack.pop()
        if h == len(strides):
            terminals += 1
            y = pow(q, n, 11**13)
            parts = (slopes[0]*y, slopes[1]*y, slopes[2]*y+1)
            assert any(sum(x % (11**j) for x in parts) >= 11**j for j in range(1, 14))
            continue
        bit = int(bool(z & (binary << h)))
        for digit, child in transitions[state].items():
            if digit == bit:
                stack.append((h+1, child, n, z))
            elif n+strides[h] <= bound:
                stack.append((h+1, child, n+strides[h], z*jumps[h] % high))
    print(q, "exponents rejected by at most thirteen base-11 digits", terminals)
    return terminals


if __name__ == "__main__":
    cases = small_bases()
    assert [lift(q, *hit) for (_, q), hit in sorted(cases.items())] == [
        8935, 13436, 8854, 8880]
