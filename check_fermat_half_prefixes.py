"""Exhaust the candidate third parts and their binary extensions at exponent 4096."""

import sys
from collections import Counter, defaultdict

P = 65537
K = 14
X = P - 2**K


def seed_prefixes(b):
    n = b + K
    modulus = 1 << n
    mask = modulus-1
    width = max(15, (n+1)//2)
    # T = 1 + 2^14 mod 2^15. Read its low digits, discarding any binomial carry.
    prefixes = [(1 << K) + 1]
    for position in range(15, width):
        place = 1 << position
        prefixes = [v for t in prefixes for v in (t, t + place)
                    if not (X*v & (v << K) & place)]
    block = 1 << width
    odd_modulus = 13*P
    inv = pow(block, -1, odd_modulus)
    limit = (modulus-1)//P
    step = block*odd_modulus
    classes = [13*(r*pow(2, b-2, P)*pow(13, -1, P) % P) for r in range(5)]
    for prefix in prefixes:
        z0 = pow(P*prefix, -1, block)
        t0 = pow(P*z0, -1, modulus)
        # If z = z0 + block*i, T changes linearly in i since block^2 = 0 mod 2^n.
        increment = -P*block*t0*t0 & mask
        stride = increment*odd_modulus & mask
        for residue in classes:
            i = (residue-z0)*inv % odd_modulus
            t = (t0 + i*increment) & mask
            for z in range(z0 + block*i, limit+1, step):
                if not (X*t & (t << K) & mask):
                    yield z, t
                t = (t + stride) & mask


def extinguish(b, z, prefix):
    scale = 1 << (b+K)
    live = [prefix]
    for position in range(40):
        if not live:
            return position
        bit = 1 << position
        high = bit*scale
        following = []
        for low in live:
            for t in (low, low+high):
                if X*t & (t << K) & high:
                    continue
                if ((P*z*t-1)//scale) & bit and not t & bit:
                    continue
                following.append(t)
        live = following
    assert not live, ("binary extension survived", b, z)
    return 40


def check_by_valuations(b, z):
    # Choose digits of B=(P*z*T-1)/S instead. The full row and its cut must
    # both have zero Legendre defect at the modulus under consideration.
    R = 1 << b
    S = R << K
    low_B = [0]
    for i in range(40):
        modulus = S << (i+1)
        inv = pow(P*z, -1, modulus)
        following = []
        for old in low_B:
            for B in (old, old + (1 << i)):
                t = (1+S*B)*inv % modulus
                parts = (X*R*t, P*z*t, (S-P*z)*t+1)
                split = (X*t, (1 << K)*t)
                if all(sum(row)//modulus == sum(v//modulus for v in row)
                       for row in (parts, split)):
                    following.append(B)
        if not following:
            return i+1
        low_B = following
    raise AssertionError(("valuation layers left a completion", b, z))


EXPECTED = {b: (0, 0) for b in range(4, 53, 4)}
EXPECTED.update({56: (2, 5), 60: (16, 8), 64: (74, 40), 68: (392, 26)})


def check_fermat_quotients():
    modulus = P*P
    assert (pow(2, P-1, modulus)-1)//P == P-4096
    for b in range(0, 96, 4):
        for digit in range(5):
            # Choose a unit T satisfying the first cut at R*T/4 = 1+P*(digit-3).
            r = pow(2, b, modulus)
            w = 4*(1 + P*(digit-3)) % modulus
            t = w*pow(r, -1, modulus) % modulus
            assert (X*w-1)//P % P == digit
            qt = (pow(t, P-1, modulus)-1)//P
            # Q(T) = 4096*j*Q(q) when T=q^(4096*j).
            jq = qt*pow(4096, -1, P) % P
            assert (b-2-jq + 16*(digit-3)) % P == 0


def main():
    check_fermat_quotients()
    compare = sys.argv[1:] == ["--compare-seeds"]
    assert compare or len(sys.argv) == 1
    reference = defaultdict(set)
    if compare:
        for line in sys.stdin:
            b, z = map(int, line.split())
            assert b in EXPECTED
            reference[b].add(z)
    for b, expected in EXPECTED.items():
        seeds = dict(seed_prefixes(b))
        if compare:
            assert set(seeds) == reference[b], b
        depths = Counter()
        for z, t in seeds.items():
            depth = extinguish(b, z, t)
            assert check_by_valuations(b, z) == depth
            depths[depth] += 1
        result = (len(seeds), max(depths, default=0))
        assert result == expected, (b, result)
        print(f"b={b}: {result[0]} candidate parts, ruled out by {result[1]} further bits", flush=True)


if __name__ == "__main__":
    main()
