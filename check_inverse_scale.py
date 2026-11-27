"""Independent integer controls for inverse exponents and the six residual families."""

from math import comb, factorial, gcd, isqrt


def primes(bound):
    return [q for q in range(2, bound + 1)
            if all(q % t for t in range(2, isqrt(q) + 1))]


def valuation(n, q):
    assert n > 0
    e = 0
    while n % q == 0:
        n //= q
        e += 1
    return e


def multiplicative_order(unit, modulus):
    assert gcd(unit, modulus) == 1
    e, residue = 1, unit % modulus
    while residue != 1:
        e += 1
        residue = residue * unit % modulus
    return e


def factorial_layers():
    layers = 0
    for B in (3, 5, 8, 12, 17, 25):
        for J in (0, 2):
            period = factorial((J + 3) * B * B)
            for p in primes(B):
                u, b = 1, p
                while b <= B:
                    for q in primes(B):
                        if q == p:
                            continue
                        E = 0
                        while q ** (E + 1) <= 2 * B:
                            E += 1
                        length = b - b // p
                        modulus = q ** (length * (J + E + 1))
                        phi = modulus - modulus // q
                        assert period % phi == 0
                        # Reduce before powering. No gigantic weight is constructed.
                        assert pow(p, (period - u) % phi, modulus) == pow(b, -1, modulus)
                        for d in range(2, b + 1):
                            n_low = (1 + d * pow(b, -1, modulus)) % modulus
                            assert n_low and valuation(n_low, q) == valuation(b + d, q)
                        for level in range(1, J + E + 2):
                            w = q ** (length * level)
                            scale = pow(b, -1, w)
                            for x in range(1, b):
                                # At fixed x, the largest admissible y minimizes its residue.
                                y = b - x
                                assert (x * scale) % w + (y * scale) % w > w
                                assert (x * scale) // w + (y * scale) // w < (b * scale) // w
                                layers += 1
                    u, b = u + 1, b * p
    return layers


def exact_multinomials():
    # Test the needed congruences just at possible common slope primes. This
    # also admits small exponents where the binomial products fit in memory.
    rows = divisible = 0
    for p, u in ((2, 1), (2, 2), (2, 3), (2, 4), (2, 5), (3, 1), (3, 2), (3, 3), (5, 1), (5, 2), (7, 1)):
        b = p ** u
        for d in range(2, b + 1):
            # For J=0 only primes dividing b+d and a possible common slope need
            # a layer. Enforce p^(e+u) = 1 at other possible slope primes as well.
            moduli = []
            for q in primes((d // 2)):
                if q == p:
                    continue
                modulus = q
                if (b + d) % q == 0:
                    layers = valuation(b + d, q)
                    modulus = q ** (layers * multiplicative_order(q, b))
                    # One more q digit than v_q(b+d) fixes the valuation of N.
                    modulus = max(modulus, q ** (layers + 1))
                moduli.append(modulus)
            for e in range(1, 11):
                if d * p ** e > 1000 or any(pow(p, e + u, w) != 1 for w in moduli):
                    continue
                Q, n = p ** e, 1 + d * p ** e
                for x in range(1, d):
                    for y in range(x, d - x + 1):
                        a, bpart = Q * x, Q * y
                        coeff = comb(n, a) * comb(n - a, bpart)
                        rows += 1
                        if coeff % p:
                            assert coeff % n == 0
                            divisible += 1
    return rows, divisible


def gap_controls():
    cases = [(5, 3, 1, 0, 0), (7, 2, 1, 0, 2), (7, 2, 1, 1, 0),
             (7, 2, 2, 0, 1), (7, 3, 1, 0, 1), (7, 5, 1, 0, 0)]
    tested = excluded = progressions = dominant = 0
    for offset, p, gap, b, v in cases:
        k, s, t = p ** gap, p ** b, p ** v
        L = (k + 1) * s
        for H in range(1, 19):
            X, G = p ** H, p ** H + 1
            D = L * G + t
            assert D == (k + 1) * p ** (H + b) + offset
            rest, factors = G, []
            q = 2
            while q * q <= rest:
                if rest % q == 0:
                    factors.append(q)
                    rest //= q ** valuation(rest, q)
                q = 3 if q == 2 else q + 2
            if rest > 1:
                factors.append(rest)
            R = 1
            for q in factors:
                if (k + 1) % q:
                    R *= q
            m0 = (L * pow(t, -1, R)) % R if R > 1 else 0
            if m0 % 2 == 0:
                m0 += R
            # Avoid expanding exponents proportional to G.
            assert m0 % 2 and (t * m0 - L) % R == 0
            for m in (*range(1, 64, 2), m0, m0 + 2 * R):
                if m * H <= v:
                    continue
                # Modulo G^2 suffices for all the valuation obstructions in (2).
                scale = pow(p, m * H - v, G * G)
                n_low = (1 + D * scale) % (G * G)
                assert n_low % G == 0
                assert (t * (n_low // G) - t * m + L) % G == 0
                for q in factors:
                    if (t * m - L) % q == 0:
                        power = q ** valuation(G, q)
                        assert n_low % (power * q) == 0
                        assert (s * G * scale) % (power * q) != 0
                        excluded += 1
                if (p, k, s, t) == (2, 2, 1, 4) and H >= 2:
                    for q in factors:
                        if q ** 3 <= G ** 2:
                            continue
                        g = G // q
                        z = (q - g) // 4
                        residue_parts = [((part % (q*q)) // q) for part in
                                         (2 * G * scale, G * scale, 4 * scale + 1)]
                        assert residue_parts == [2*z, z, (m*g) % q]
                        assert (sum(residue_parts) >= q) == ((m*g) % q >= q - 3*z)
                        dominant += 1
                if (t * m - L) % R == 0:
                    assert all((k + 1) % q == 0 or n_low % (G * q) == 0
                               for q in factors)
                    progressions += 1
                tested += 1
    return tested, excluded, progressions, dominant


def divisor_and_large_controls():
    bounded = huge = 0
    for G in range(2, 12001):
        # Build the interval of candidate divisors independently of prime grouping.
        found = any(G % z == 0 and gcd(z, G // z) == 1 and z ** 3 >= G
                    for z in range(2, isqrt(G) + 1))
        if not found:
            n = G
            atoms = []
            for q in primes(isqrt(G)):
                if n % q == 0:
                    power = q ** valuation(n, q)
                    atoms.append(power)
                    n //= power
            if n > 1:
                atoms.append(n)
            assert max(atoms) ** 3 > G ** 2
        bounded += 1
    for p, k, s, t in [(3, 3, 1, 1), (2, 2, 1, 4), (2, 2, 2, 1),
                       (2, 4, 1, 2), (3, 3, 1, 3), (5, 5, 1, 1)]:
        v, L = valuation(t, p), (k + 1) * s
        for H in (257, 512, 1001):
            G = p ** H + 1
            R = G
            while gcd(R, k + 1) > 1:
                R //= gcd(R, k + 1)
            m = (L * pow(t, -1, R)) % R if R > 1 else 0
            if m % 2 == 0:
                m += R
            Q = pow(p, m * H - v, G * G)
            assert (1 + (L * G + t) * Q) % (G * R) == 0
            if (p, k, s, t) == (2, 2, 1, 4):
                m = 2 ** (H - 2) + 1
                Q = pow(2, m * H - 2, G * G)
                assert (1 + (3 * G + 4) * Q) % (G * G) == 0
            huge += 1
    return bounded, huge


if __name__ == "__main__":
    print("Inverse residue floor comparisons:", factorial_layers())
    print("Exact multinomials, including p-free rows:", exact_multinomials())
    print("Residual exponents, excluded factors, class hits, dominant-prime layers:", gap_controls())
    print("Unitary divisor checks, large inner-exponent certificates:", divisor_and_large_controls())
