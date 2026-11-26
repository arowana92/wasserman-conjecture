"""Check the two-prime gap criterion directly against exact coefficients."""

from itertools import combinations
from math import comb, gcd, isqrt


def prime_powers(n):
    powers = []
    m = n
    for p in range(2, isqrt(n) + 1):
        if m % p:
            continue
        w = 1
        while m % p == 0:
            m //= p
            w *= p
            powers.append((p, w))
    if m > 1:
        powers.append((m, m))
    return powers


def top_power(n, p):
    w = p
    while w*p <= n:
        w *= p
    return w


def block_fails(n, p_power, q, q_power, gap):
    # A shared avoiding row would contain two multiples of p_power*q_power.
    units = gap // (p_power * q_power)
    w = q
    while w <= n:
        s = (-p_power) % w
        if s*units < 2*w - (n//q_power) % w:
            return True
        w *= q
    return False


def check_coefficients(limit):
    obstructed = joint_avoiders = split_tests = 0
    for n in range(3, limit + 1):
        pairs = [(p, P, q, Q) for (p, P), (q, Q) in combinations(prime_powers(n), 2)
                 if p != q and 2*min(top_power(n, p), top_power(n, q)) > n]
        if not pairs:
            continue
        rows = [((a, b, n-a-b), comb(n, a)*comb(n-a, b))
                for a in range(1, n//3+1) for b in range(a, (n-a)//2+1)]
        for p, P, q, Q in pairs:
            left = n-top_power(n, p)
            right = n-top_power(n, q)
            singles_p = [(row, k) for row, k in rows if k % p]
            singles_q = [(row, k) for row, k in rows if k % q]
            assert all(sum(row[:2]) <= left for row, _ in singles_p)
            assert all(sum(row[:2]) <= right for row, _ in singles_q)
            fail = (block_fails(n, P, q, Q, min(left, right)) or
                    block_fails(n, Q, p, P, min(left, right)))
            same = any(k % q for _, k in singles_p)
            if same:
                joint_avoiders += 1
            if fail:
                assert not same, (n, P, Q)
                obstructed += 1
            if (left*right)**2 < 16*n*n*(n-1):
                for _, a in singles_p:
                    for _, b in singles_q:
                        assert gcd(n*(n-1), gcd(a, b)) > 1
                        split_tests += 1
    return obstructed, joint_avoiders, split_tests


def check_mixed_blocks(limit):
    # Inspect actual entries and their next blocks for all allowed pairs of low residues.
    checked = rejected = compatible = 0
    powers = sorted((w, p) for w in range(2, limit+1)
                    for p, v in prime_powers(w) if v == w)
    for (P, p), (Q, q) in combinations(powers, 2):
        s = Q-P
        if s < 2 or p == q:
            continue
        ip, iq = pow(s, -1, P), pow(s, -1, Q)
        threshold = min(ip, P-ip, iq, Q-iq)
        modulus = P*Q
        # One CRT representative for each allowed residue pair, then successive periods.
        parts = []
        for delta in (0, 1):
            for eta in (0, 1):
                first = delta + P*((eta-delta)*pow(P, -1, Q) % Q)
                for a in range(first, modulus*(threshold+2)//s + 2*modulus, modulus):
                    if a >= 2:
                        parts.append((a, delta, eta))
        for i, (a, da, ea) in enumerate(parts):
            for b, db, eb in parts[i:]:
                if da+db > 1 or ea+eb > 1:
                    continue
                checked += 1
                fits = ((a//P)%P + (b//P)%P < P and (a//Q)%Q + (b//Q)%Q < Q)
                if s*(a+b)+P < modulus*threshold:
                    assert not fits, (P, Q, a, b)
                    rejected += 1
                if fits:
                    compatible += 1
    return checked, rejected, compatible


def mixed_interval(P, Q):
    s = Q-P
    step = P*Q
    low, high = Q**4, Q**4 + Q**3//(4*s*s)
    counts = []
    for delta, eta in ((0, 0), (0, 1), (1, 0), (1, 1)):
        residue = delta + P*((eta-delta)*pow(P, -1, Q) % Q)
        first = low + (residue-low) % step
        last = high - (high-residue) % step
        if first > last:
            counts.append(0)
            continue
        assert first % P == last % P == delta and first % Q == last % Q == eta
        h, k = last-P**4, last-Q**4
        assert (h*k)**2 < 16*last**2*(last-1)
        ip, iq = pow(s, -1, P), pow(s, -1, Q)
        assert s*k+P < step*min(ip, P-ip, iq, Q-iq)
        counts.append(1+(last-first)//step)
    return counts


def check_consecutive(P, Q):
    step = P*Q
    q = prime_powers(Q)[0][0]
    assert P+1 == Q and Q >= 10
    top = Q**3//4
    special = set()
    for delta in (0, 1):
        for eta in (0, 1):
            for z in range(2+top//step):
                part = step*z+Q*delta-P*eta
                if not 0 < part <= top:
                    continue
                u = z-delta+eta
                assert (part//Q)%Q == (-u)%Q
                if u <= 0:
                    special.add(part)
    assert special == {1, Q, Q**2}
    for m in range(1, top//step+1):
        n = Q**4+step*m
        cofactor = m
        while cofactor % q == 0:
            cofactor //= q
        v = m//cofactor
        assert n%(Q*v) == 0 and (n//(Q*v))%q
        S = gcd(n-1, Q**2-1)
        assert S == P*gcd(m, Q+1) and S <= P*cofactor
        assert (Q*v)**3 * S**2 < 27*(n-1)**2
    return top//step


def check_cubic(P, Q):
    s = Q-P
    assert s >= 2 and (4*s)**4 <= Q
    high = Q**3+isqrt(isqrt(Q**9))
    h, k = high-P**3, high-Q**3
    # Compare the largest left sides to the right sides at the lower endpoint.
    assert 2*P**3 > high and (h*k)**2 < 16*Q**6*(Q**3-1)
    ip, iq = pow(s, -1, P), pow(s, -1, Q)
    assert (s*k+P) < P*Q*min(ip, P-ip, iq, Q-iq)
    return k


def quartics(P, Q):
    # The upper half of the interval supplies examples at its outer scale.
    s = Q-P
    assert gcd(P, Q) == 1 and 1 <= 10*s <= Q
    step = P*Q
    first = (Q**4 + step-1)//step
    last = (Q**4 + Q**3//(4*s))//step
    n = last*step
    assert first <= last and 2*P**4 > n
    h, k = n-P**4, n-Q**4
    assert (h*k)**2 < 16*n*n*(n-1)
    assert s*(k//step) <= Q
    return first*step, n, last-first+1


if __name__ == '__main__':
    print('Small-weight checks (obstructions, pairs of powers with a joint avoider, '
          'split-row comparisons):', check_coefficients(510))
    for P, Q in ((29, 31), (1009, 1013), (2**19, 3**12), (1000000007, 1000000009)):
        print(P, Q, *quartics(P, Q), mixed_interval(P, Q))
    print('Mixed blocks (tested, rejected by the bound, compatible):', check_mixed_blocks(100))
    print('Consecutive-power multipliers:', [check_consecutive(P, P+1)
          for P in (16, 31, 256, 65536)])
    print('Cubic interval width:', check_cubic(1000000007, 1000000009))
