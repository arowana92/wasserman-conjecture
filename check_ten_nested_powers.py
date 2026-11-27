"""Local carry certificates when Q is one modulo P at cofactor ten."""

from itertools import product
from math import gcd

from check_static_power_transfer import prime_powers
from check_three_factors import valuation


def shapes(length):
    for u, x in ((3, 3), (3, 6), (3, 9), (7, 7)):
        for k in range(1, 10-x):
            t = 10-x-k
            if gcd(k, 5*u) != 1:
                continue
            P = -pow(k, -1, 5*u) % (5*u)

            def extend(digits, previous):
                if len(digits) == length:
                    m = sum(z*P**i for i, z in enumerate(digits))
                    Q = 1+m*P
                    if (sum(digits) % 2 == 0 and (1+10*P*Q) % u == 0
                            and x*P*Q % 5 == 1):
                        yield u, x, k, t, digits
                    return
                for z in range(1, t*previous+1):
                    yield from extend(digits+(z,), z)

            yield from extend((), 1)


def entries(shape, P):
    u, x, k, t, digits = shape
    m = sum(z*P**i for i, z in enumerate(digits))
    Q = m*P+1
    n = 10*P*Q+1
    return Q, n, (x*P*Q, P*(t*Q-m), Q*(k*P+1))


def allowed(shape, P, p, modulus):
    u, x, k, t, _ = shape
    Q, n, row = entries(shape, P)
    free = sum(a % modulus for a in row) < modulus
    if p == 2:
        if P*Q % 2 == 0:
            return False
        if u == 3:
            return free
        v = n * pow(7, -1, modulus) % modulus
        return free or sum((a*v) % modulus for a in (1, 2, 4)) < modulus
    if (k*P+1) % p or not free:
        return False
    if p == 5:
        return x*P*Q % 5 == 1
    return n % u == 0 and (modulus < u*u or n % (u*u) != 0)


def first_empty(shape, p, depth, cap):
    good = [0]
    modulus = 1
    for pos in range(1, depth+1):
        prev = modulus
        modulus *= p
        good = [z for base in good for z in range(base, modulus, prev)
                if allowed(shape, z, p, modulus)]
        if not good:
            return pos
        if len(good) > cap:
            break
    return None


def independent_shapes(length):
    """Try every digit box and all residue classes for P, without inverting k."""
    found = set()
    for x in range(3, 10):
        for u in (3, 7):
            if x % u:
                continue
            for t in range(1, 10-x):
                k = 10-x-t
                roots = [a for a in range(1, 10*u, 2) if (k*a+1) % (5*u) == 0]
                for digits in product(*(range(1, t**j+1) for j in range(1, length+1))):
                    if any(b > t*a for a, b in zip(digits, digits[1:])):
                        continue
                    for a in roots:
                        m = 0
                        for z in reversed(digits):
                            m = a*m+z
                        W = a*(1+a*m)
                        if W % 2 and x*W % 5 == 1 and (10*W+1) % u == 0:
                            found.add((u, x, k, t, digits))
    return found


def check_obstruction(shape, p, depth):
    """Rebuild a claimed empty tree from individual base-p digits of exact parts."""
    u, x, k, t, digits = shape
    assert u == 3  # All candidates in the three checked boxes have this scale.
    roots = [0]
    place = 1
    for j in range(1, depth+1):
        next_roots = []
        for digit in range(p):
            for base in roots:
                P = base + digit*place
                m = 0
                for z in reversed(digits):
                    m = P*m+z
                Q = 1 + m*P
                n = 1 + 10*P*Q
                if p == 2 and P*Q % 2 == 0:
                    continue
                if p == 3 and ((P*k+1) % 3 or n % 3 or (j >= 2 and n % 9 == 0)):
                    continue
                if p == 5 and ((P*k+1) % 5 or x*P*Q % 5 != 1):
                    continue
                a = x*P*Q
                b = Q*(1+k*P)
                c = n-a-b
                if sum((v//place) % p for v in (a, b, c)) < p:
                    next_roots.append(P)
        roots = next_roots
        place *= p
    assert not roots, (shape, p, depth)


def large_powers():
    special = (3, 6, 2, 2, (2, 4, 4))
    for length in range(1, 4):
        candidates = set(shapes(length))
        assert candidates == independent_shapes(length)
        assert len(candidates) == (0, 9, 470)[length-1]
        maxima = {}
        for shape in sorted(candidates):
            if shape == special:
                continue
            for p in (2, shape[0], 5):
                depth = first_empty(shape, p, 60, 1200)
                if depth:
                    check_obstruction(shape, p, depth)
                    maxima[p] = max(maxima.get(p, 0), depth)
                    break
            else:
                raise AssertionError(shape)
        print(f"{length} digits: {len(candidates)} shapes, last empty layers {maxima}")
    for P in range(1, 5000):
        Q, n, (_, b, c) = entries(special, P)
        T = 2*P+1
        assert b == 4*P**3*T and b+c == T*T*(4*P*P+1)
        if T % 5 == 0:
            assert valuation(n, (n-b-c, b, c), 5) > 0


def small_powers():
    completed = [0]*3
    largest_position = -1
    for P, p in prime_powers(500).items():
        if p in (2, 5):
            continue
        for u, x in ((3, 3), (3, 6), (7, 7)):
            for k in range(1, 10-x):
                t = 10-x-k
                if (k*P+1) % (5*u):
                    continue
                # A next occupied position past the cutoff gives opposite nonzero
                # digits in the final two upper parts, regardless of later digits.
                stack = [(0, 1, 0, 0)]
                while stack:
                    low, place, used, pos = stack.pop()
                    if place > 10*(1+P*low):
                        continue
                    modulus = place*P
                    for digit in range(P):
                        m = low + digit*place
                        upper = (x+x*P*m, t+(t*P-1)*m, k+(k*P+1)*m)
                        if sum(v % modulus for v in upper) >= modulus:
                            continue
                        count = used + bool(digit)
                        if count < 3:
                            stack.append((m, modulus, count, pos+1))
                        if digit == 0:
                            continue
                        completed[count-1] += 1
                        largest_position = max(largest_position, pos)
                        Q = 1+m*P
                        W = P*Q
                        n = 1+10*W
                        if not Q % 2 or x*W % 5 != 1 or n % u or n % (u*u) == 0:
                            continue
                        row = (x*W, P*(t*Q-m), Q*(k*P+1))
                        if any(valuation(n, row, r) for r in (5, u, p)):
                            continue
                        assert valuation(n, row, 2) > 0
                        assert u == 3 or valuation(n, tuple(a*(n//7) for a in (1, 2, 4)), 2) > 0
    print(f"P <= 500: {completed} sparse quotients, last occupied position "
          f"{largest_position}, every remaining row excluded by valuations")


def independent_small_powers():
    """Place the nonzero digits using the gap theorem and budget each part separately."""
    def free(parts, prime):
        while any(parts):
            if sum(a % prime for a in parts) >= prime:
                return False
            parts = [a//prime for a in parts]
        return True

    checked = 0
    for P, p in prime_powers(500).items():
        if p in (2, 5):
            continue
        for u, x in ((3, 3), (3, 6), (7, 7)):
            for t in range(1, 10-x):
                k = 10-x-t
                if (1+k*P) % (5*u):
                    continue
                assert P >= 10-x  # The base-p digit-gap bound is at most 2 log_p(P).
                stack = [(0, -1, 0)]
                while stack:
                    low, prev, used = stack.pop()
                    for pos in (prev+1, prev+2):
                        place = P**pos
                        Q0 = 1+P*low
                        # At this position, the last two upper parts change by -z,+z.
                        target = (10*Q0//place) % P
                        first = (x*Q0//place) % P
                        second = ((t*Q0-low)//place) % P
                        for b in range(target-first+1):
                            z = (second-b) % P
                            if not z:
                                continue
                            m = low+z*place
                            if used < 2:
                                stack.append((m, pos, used+1))
                            Q = 1+P*m
                            weight = 1+10*P*Q
                            if Q % 2 == 0 or weight % u or weight % (u*u) == 0:
                                continue
                            a, c = x*P*Q, Q*(k*P+1)
                            if a % 5 != 1:
                                continue
                            checked += 1
                            row = (a, weight-a-c, c)
                            if all(free(row, r) for r in (p, u, 5)):
                                assert not free(row, 2)
                                assert u == 3 or not free([i*(weight//7) for i in (1, 2, 4)], 2)
    print(f"Independent sparse-position enumeration excluded {checked} static candidates")


def binary_transitions():
    slopes = (6, 3, 1)
    table = {
        (3, 1, 1): (None, (4, 2, 1)),
        (4, 2, 1): ((2, 1, 0), (5, 2, 1)),
        (2, 1, 0): ((1, 0, 0), (4, 2, 0)),
        (5, 2, 1): (None, None),
        (1, 0, 0): ((0, 0, 0), None),
        (4, 2, 0): ((2, 1, 0), None),
        (0, 0, 0): ((0, 0, 0), None),
    }
    for state, choices in {**table, (0, -1, 1): (None, (3, 1, 1))}.items():
        for bit in (0, 1):
            if state == (0, -1, 1) and bit == 0:
                continue  # P is odd.
            values = tuple(a*bit+c for a, c in zip(slopes, state))
            actual = tuple(v//2 for v in values) if sum(v%2 for v in values) <= 1 else None
            assert actual == choices[bit]
    print("binary automaton for the even Q case verified")


if __name__ == "__main__":
    binary_transitions()
    large_powers()
    small_powers()
    independent_small_powers()
