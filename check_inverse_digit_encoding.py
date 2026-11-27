"""Check the delayed-digit injection and the sieve for three predecessor primes."""
from itertools import combinations_with_replacement
from math import comb, gcd, prod


def children(low, place, earlier_place, p, W, s, k):
    """Criterion on the next digit of T. Here place = W*earlier_place."""
    parent = low//earlier_place % p
    inv = pow(s, -1, p*earlier_place)
    result = []
    for digit in range(p):
        value = low+digit*place
        numerator = k*value-s
        assert numerator % W == 0
        x = (numerator//W)*inv % (p*earlier_place)
        if x//earlier_place <= parent:
            result.append(digit)
    assert len(result) == parent+1
    return result


def encode(t, p, e, length, s, k):
    """Return a word satisfying the distance-e inequalities, or reject."""
    W = p**e
    word = []
    place = 1
    low = 0
    for j in range(length):
        digit = t//place % p
        if j < e:
            out = digit
        else:
            allowed = children(low, place, place//W, p, W, s, k)
            if digit not in allowed:
                return None
            out = p-len(allowed)+allowed.index(digit)
            assert out >= digit and out+word[j-e] >= p-1
        word.append(out)
        low += digit*place
        place *= p
    return tuple(word)


def decode(word, p, e, s, lift):
    W = p**e
    t = sum(v*p**j for j, v in enumerate(word[:e]))
    if len(word) <= e:
        return t
    k = s*pow(t, -1, W) % W + lift*W
    place = W
    for value in word[e:]:
        allowed = children(t, place, place//W, p, W, s, k)
        rank = value+len(allowed)-p
        assert 0 <= rank < len(allowed)
        t += place*allowed[rank]
        place *= p
    return t


def digit_le(x, n, p):
    while x:
        if x % p > n % p:
            return False
        x //= p
        n //= p
    return True


def test_encoding():
    tries = accepted = 0
    for p, length in ((2, 12), (3, 7), (5, 5), (7, 4)):
        for e in range(1, length):
            W, bound = p**e, p**length
            for s in (1, 2, 3, 4, 5):
                if s % p == 0:
                    continue
                for t in range(1, bound):
                    if t % p == 0:
                        continue
                    first = s*pow(t, -1, W) % W
                    for lift in range(s):
                        k = first+lift*W
                        actual = (k*t-s)//W * pow(s, -1, bound//W) % (bound//W)
                        expect = digit_le(actual, t, p)
                        word = encode(t, p, e, length, s, k)
                        assert (word is not None) == expect
                        tries += 1
                        if expect:
                            assert decode(word, p, e, s, lift) == t
                            accepted += 1
    print("prefix tests", tries, "accepted and recovered", accepted)


def prefix_tree(p, exponent):
    # Classes N = 1 mod p. The class 1 modulo p^exponent has unbounded valuation.
    counts = [1]
    for e in range(1, exponent):
        W = p**e
        length = exponent-e
        bound = p**length
        if W <= 4:
            counts.append(0)
            continue
        if length < e:
            counts.append((p-1)*(bound//p))
            continue
        def extend(low, place, k):
            if place == bound:
                return 1
            return sum(extend(low+d*place, p*place, k)
                       for d in children(low, place, place//W, p, W, 1, k))
        count = 0
        for low in range(1, W):
            if low % p:
                k = pow(low, -1, W)
                if k >= 3:
                    count += extend(low, W, k)
        counts.append(count)
    return counts


def direct_sieve(p, exponent):
    modulus = p**exponent
    counts = [1]+[0]*(exponent-1)
    for n in range(1+p, modulus, p):
        t = n-1
        e = 0
        while t % p == 0:
            t //= p
            e += 1
        W = p**e
        if W <= 4:
            continue
        known = modulus//W
        if known < W:
            counts[e] += 1
            continue
        k = pow(t, -1, W)
        if k >= 3 and digit_le(k*t % known, n % known, p):
            counts[e] += 1
    return counts


def chain_bound(p, exponent):
    """Upper bound from the independent encoded chains, including the first digits."""
    paths = [[0]*p, [1]*p]
    for length in range(2, exponent+1):
        paths.append([sum(paths[length-1][p-1-first:]) for first in range(p)])
    counts = [1]
    for e in range(1, exponent):
        W = p**e
        length = exponent-e
        if W <= 4:
            counts.append(0)
            continue
        if length < e:
            counts.append((p-1)*p**(length-1))
            continue
        chains = [paths[1+(length-1-i)//e] for i in range(e)]
        total = sum(chains[0][1:])*prod(sum(tail) for tail in chains[1:])
        for k in (1, 2):
            if k % p == 0:
                continue
            initial = pow(k, -1, W)
            for_start = 1
            for options in chains:
                for_start *= options[initial % p]
                initial //= p
            total -= for_start
        counts.append(total)
    return counts


def check_chain_counts():
    for p, exponent in ((2, 12), (3, 8), (5, 6)):
        predicted = chain_bound(p, exponent)
        for e in range(1, exponent):
            W = p**e
            length = exponent-e
            accepted = 0
            for word in range(p**length):
                if W <= 4 or word % p == 0:
                    continue
                if length >= e and pow(word, -1, W) < 3:
                    continue
                digits = [(word//p**i) % p for i in range(length)]
                accepted += all(a+b >= p-1 for a, b in zip(digits, digits[e:]))
            assert accepted == predicted[e]


def binary_fibonacci_bound(exponent):
    fibonacci = [0, 1]
    for _ in range(exponent+2):
        fibonacci.append(fibonacci[-1]+fibonacci[-2])
    total = 1
    for e in range(3, exponent):
        length = exponent-e
        if length < e:
            total += 2**(length-1)
            continue
        short, longer = divmod(length, e)
        # The first chain starts at 1. The other initial digits cannot all be 0.
        first_length = short+bool(longer)
        other_long = max(0, longer-1)
        other_short = e-1-other_long
        choices = fibonacci[short+3]**other_long * fibonacci[short+2]**other_short
        all_zero = fibonacci[short+1]**other_long * fibonacci[short]**other_short
        total += fibonacci[first_length+1]*(choices-all_zero)
    return total


def valuation(value, p):
    exponent = 0
    while value % p == 0:
        value //= p
        exponent += 1
    return exponent


def progressions():
    comparisons = 0
    for p in (2, 3, 5):
        for s in (1, 2, 3):
            if s % p == 0:
                continue
            for step in (1, 2, 3, 6, p**3, 2*p**5, 7*p**8):
                f = valuation(step, p)
                for first in (4, 1+p, 1+p**2, 1+s*p**4, 1+s*step):
                    weights = [first+i*step for i in range(79) if (first+i*step-1) % s == 0]
                    if not weights:
                        continue
                    n = 0
                    while p**n < len(weights):
                        n += 1
                    a = min(f, valuation(first-1, p))
                    tail = 0
                    groups = {}
                    for weight in weights:
                        e = valuation(weight-1, p)
                        if not e:
                            continue
                        if a == f and e-f >= n:
                            tail += 1
                            continue
                        frozen = f-a if a < f else 0
                        length = n if a < f else n-(e-f)
                        W = p**e
                        t = (weight-1)//W
                        block = t//p**frozen % p**length
                        blocks, words = groups.setdefault(e, (set(), set()))
                        assert block not in blocks
                        blocks.add(block)
                        for lift in range(s):
                            k = s*pow(t, -1, W) % W + lift*W
                            part = k*(t//s)
                            if not digit_le(part, weight, p):
                                continue
                            word = []
                            for j in range(length):
                                place = p**(frozen+j)
                                digit = t//place % p
                                if j >= e:
                                    allowed = children(t % place, place, place//W, p, W, s, k)
                                    digit = p-len(allowed)+allowed.index(digit)
                                    assert digit+word[j-e] >= p-1
                                word.append(digit)
                            key = (lift, tuple(word))
                            assert key not in words
                            words.add(key)
                            comparisons += 1
                    assert tail <= 1
    print("Progression block injections", comparisons)


def exact_rows():
    checks = 0
    for p in (2, 3, 5, 7):
        for n in range(1+p, 130, p):
            t = n-1
            e = 0
            while t % p == 0:
                t //= p
                e += 1
            W = p**e
            for a, b in combinations_with_replacement(range(1, n//2+1), 2):
                c = n-a-b
                if c < b or comb(n, a)*comb(n-a, b) % p == 0:
                    continue
                shifted, = [x for x in (a,b,c) if x % W == 1]
                s = t//gcd(shifted, t)
                k = shifted*s//t
                assert 1 <= k < s*W and (k*t-s) % W == 0
                word = encode(t, p, e, 10, s, k)
                assert word is not None and decode(word, p, e, s, k//W) == t % p**10
                checks += 1
    print("exact trinomial comparisons", checks)


if __name__ == "__main__":
    test_encoding()
    check_chain_counts()
    progressions()
    exact_rows()
    for p, exponent in ((2,24), (3,13), (5,9)):
        counts = prefix_tree(p, exponent)
        assert counts == direct_sieve(p, exponent)
        assert all(a <= b for a, b in zip(counts, chain_bound(p, exponent)))
        print("modulus", p, exponent, "surviving classes", sum(counts),
              "of", p**(exponent-1), "by valuation", counts, flush=True)
    for p, exponent in ((2,128), (3,81), (5,55)):
        limit = sum(chain_bound(p, exponent))
        if p == 2:
            assert limit == binary_fibonacci_bound(exponent)
        total = p**(exponent-1)
        print("modulus", p, exponent, "at most", limit, "of", total,
              "fraction below", (10**20*limit+total-1)//total, "/ 10^20")
