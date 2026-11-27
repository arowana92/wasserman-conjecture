"""Check the binary mask, character-cancellation and sieve identities."""
from fractions import Fraction
from math import gcd, lcm

# Positions represent (previous bit, carry of Q+2Q).
STATES = ((0,0),(1,0),(1,1),(0,1))
POTENTIAL = (0,1,2,1)
WEIGHT = (2,1,1,2)


def transition(index, bit):
    previous, carry = STATES[index]
    digit = bit + previous + carry
    if digit % 2 and previous:
        return None
    return STATES.index((bit, digit//2))


def follow(state, word, length):
    for i in range(length):
        state = transition(state, (word >> i) & 1)
        if state is None:
            return None
    return state


def binary_checks():
    for s in range(4):
        weighted = 0
        for b in (0,1):
            t = transition(s,b)
            if t is not None:
                assert 2*b-1 <= POTENTIAL[t]-POTENTIAL[s]
                weighted += (1 << b)*WEIGHT[t]
        assert weighted == 2*WEIGHT[s]
        assert any(follow(s, word, 3) == 0 for word in range(8))
    for word in (0,3,11):
        assert follow(0,word,6) == 0
    for q in range(1 << 17):
        top = q.bit_length()
        valid = follow(0,q,top+2) == 0
        assert valid == (3*q & (2*q) == 0)
        if valid and q:
            assert 4**bin(q).count("1") <= 6*q-2
    for ones in range(2,150):
        q = (2**(2*ones-1)+1)//3
        assert bin(q).count("1") == ones and 3*q & (2*q) == 0
        assert 4**ones == 6*q-2


def inverse_checks():
    # Odd composite moduli are also valid controls for the mask arithmetic.
    total = examples = 0
    for q in range(3,40000,2):
        if 3*q & (2*q):
            continue
        inv2 = (q+1)//2
        orbit = {}
        x = inv2
        while x not in orbit:
            orbit[x] = len(orbit)
            x = x*inv2 % q
        by_orbit = {k for a,k in orbit.items() if a % 15 == 0 and a & ~q == 0}
        by_mask = set()
        a = q
        while a:
            if a != q and a % 15 == 0 and gcd(a,q) == 1:
                target = pow(2*a,-1,q)
                if target in orbit:
                    k = (-1-orbit[target]) % len(orbit)
                    assert pow(2,k+1,q)*a % q == 1
                    by_mask.add(k)
            a = (a-1) & q
        assert by_mask == by_orbit and len(by_mask)**2 <= 6*q-2
        total += len(by_mask)
        examples += bool(by_mask)
    assert examples > 0
    print('Inverse residues at odd moduli below 40000:', total, 'across',examples,'moduli')


def distribution_check():
    # Compare the exact Markov masses on residue classes to all unrestricted
    # input words, including prefixes whose eventual zero padding would fail.
    for modulus in (3,11,15,33):
        counts = {(0,0):1}
        for h in range(12):
            following = {}
            for (s,r), count in counts.items():
                for b in (0,1):
                    t = transition(s,b)
                    if t is None:
                        continue
                    key = (t,(r+(b << h)) % modulus)
                    following[key] = following.get(key,0) + (count << b)
            counts = following
        expected = [0]*modulus
        for q in range(1 << 12):
            t = follow(0,q,12)
            if t is not None:
                expected[q % modulus] += (1 << bin(q).count("1"))*WEIGHT[t]
        observed = [0]*modulus
        for (s,r), count in counts.items():
            observed[r] += count*WEIGHT[s]
        assert observed == expected and sum(observed) == WEIGHT[0]*(1 << 12)


def character_checks():
    for d in range(3,180,2):
        for a in range(1,d):
            # The common shift of three digits leaves the fraction nonintegral.
            v = 3 if 3*a % d else 11
            r = v*a*8 % d
            for _ in range(d.bit_length()+1):
                if d <= 4*r <= 3*d:
                    break
                r = 2*r % d
            else:
                raise AssertionError((a,d))


def sieve_check(z):
    mu = [0]*(z+1)
    phi = [0]*(z+1)
    for n in range(1,z+1):
        mu[n] = 1 if n == 1 else -sum(mu[r] for r in range(1,n) if n%r == 0)
        phi[n] = sum(gcd(a,n) == 1 for a in range(1,n+1))
    divisors = [n for n in range(1,z+1,2) if mu[n]]
    g = sum(Fraction(1,phi[r]) for r in divisors)
    weights = {d: d*sum(Fraction(mu[r//d]*mu[r],phi[r]*g)
                        for r in divisors if not r%d) for d in divisors}
    assert weights[1] == 1 and all(abs(u) <= d for d,u in weights.items())
    main = sum(a*b/lcm(d,e) for d,a in weights.items() for e,b in weights.items())
    assert main == 1/g


if __name__ == '__main__':
    binary_checks()
    inverse_checks()
    distribution_check()
    character_checks()
    for z in (1,5,11,24,39):
        sieve_check(z)
    print('Digit, period, cancellation and sieve identities verified.')
