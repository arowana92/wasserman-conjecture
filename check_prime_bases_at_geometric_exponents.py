"""Integer checks for stable exponent digits, simultaneous cuts, and the explicit base."""

from math import comb, gcd, isqrt, prod


def factors(n):
    answer = {}
    a = 2
    while a*a <= n:
        while n % a == 0:
            answer[a] = answer.get(a, 0) + 1
            n //= a
        a += 1
    if n > 1:
        answer[n] = 1
    return answer


def prime_by_order(n):
    # Completely factor n-1 by trial division, and exhibit an element of order n-1.
    primes = factors(n - 1)
    assert prod(p**e for p, e in primes.items()) == n - 1
    assert all(all(p % a for a in range(2, isqrt(p) + 1)) for p in primes)
    for g in range(2, 100):
        if pow(g, n-1, n) == 1 and all(pow(g, (n-1)//p, n) != 1 for p in primes):
            return g, primes
    raise AssertionError("No generator among the tested candidates")


def valuation(a, p):
    v = 0
    while a % p == 0:
        a //= p
        v += 1
    return v


def low_carry(x, y, p, n):
    for _ in range(n):
        if x % p + y % p >= p:
            return True
        x //= p
        y //= p
    return False


def normalized(q, p, e, k, n):
    low = p**(k+1)
    t = pow(q, (p-1)*e*p**k, low*p**n)
    assert t % low == 1
    return t//low


def check_stability():
    layers = fixed = cuts = positive = negative = 0
    for p in (2, 3, 5, 7):
        f = 2 if p == 2 else 1
        for n in (2, 3, 4):
            h = p**n
            for offset in (7, 10**8+51):
                interval = range(offset, offset + min(h, 241))
                for e in (1, p+1, p*p+1):
                    for r in range(1, p**f):
                        if gcd(r, p) > 1:
                            continue
                        bases = [q for q in interval if q % p**f == r]
                        for j in (1, e, (p-1)*e*p**7):
                            width = p**(f + valuation(j, p))
                            words = [pow(q, j, width*h) for q in bases]
                            assert len({w % width for w in words}) <= 1
                            assert len({w//width for w in words}) == len(bases)
                            fixed += len(bases)
                        reference = [normalized(q, p, e, n, n) for q in bases]
                        assert len(set(reference)) == len(bases)
                        for q, lo in zip(bases, reference):
                            for k in (n, n+1, n+7):
                                assert normalized(q, p, e, k, n) == lo
                                layers += 1
                            for d in (3, 7, 11):
                                if d % p == 0:
                                    continue
                                k = n+7
                                shift = p**(k+1)
                                t = pow(q, (p-1)*e*p**k, shift*h)
                                for x in range(1, d):
                                    first = (x*t)//shift
                                    second = ((d-x)*t)//shift
                                    expected = low_carry(x*lo, (d-x)*lo, p, n)
                                    assert low_carry(first, second, p, n) == expected
                                    cuts += 1
                                    if expected:
                                        positive += 1
                                    else:
                                        negative += 1
    assert positive and negative
    return layers, fixed, cuts


def check_progression():
    checks = small = 0
    for p in (2, 3, 5, 7):
        for m in (2, 3, 4):
            block = p**m
            bound = block//gcd(p, 2)
            for q in (block - 1, block - 1 + block*bound):
                # The cut arithmetic also holds for composite test bases.
                for j in (1, 2, 3, 4, 6, 12, 2*p**12):
                    v = valuation(j, p)
                    shift = block*p**v
                    t = pow(q, j, shift*bound)
                    if j % 2 == 0:
                        assert t % shift == 1
                        assert t//shift == (-j//p**v) % bound
                    limit = block if j % 2 else bound//(j//p**v)
                    for d in range(2, min(limit, 32) + 1):
                        for x in range(1, d):
                            layer = block if j % 2 else shift*bound
                            assert x*t % layer + (d-x)*t % layer >= layer
                            checks += 1
                            if j <= 3 and d*q**j <= 1200:
                                assert comb(d*q**j, x*q**j) % p == 0
                                small += 1
    return checks, small


if __name__ == "__main__":
    print("stable residues, injective words, tail cut comparisons:", check_stability())
    print("progression cuts, exact binomial checks:", check_progression())
    print("generator and factorization for the explicit prime:", prime_by_order(2**61-1))
