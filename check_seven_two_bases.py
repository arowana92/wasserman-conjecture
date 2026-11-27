"""Bound both exponents of every remaining odd base pair through 47."""
from bisect import bisect_left
from collections import deque
from itertools import combinations

from flint import arb, ctx, fmpz_mat

from prime_certificates import primes_through


PASSES = ((130, 2, 82), (82, 3, 62), (62, 4, 42), (42, 5, 22), (22, 6, None))


def cyclic_closure(residues, modulus):
    queue = deque(set(residues))
    seen = set(queue)
    while queue:
        nxt = 7*queue.popleft() % modulus
        if nxt not in seen:
            seen.add(nxt)
            queue.append(nxt)
    return seen


def possible_prefixes(p, q):
    # 6(t+1) = 0 mod 72 and t = 11 mod 72. Also p*q divides t.
    # Working modulo 432 keeps the complete information from t mod 72.
    data = []
    for modulus, target in ((432, 72), (6*p*q, 6)):
        live = cyclic_closure([target], modulus)
        by_length = [None]*7
        for j in range(6, 0, -1):
            mask = bytearray(modulus)
            for n in live:
                mask[n] = 1
            by_length[j] = mask
            live = cyclic_closure(((x-1) % modulus for x in live), modulus)
        data.append((modulus, by_length))
    return data


def rounded(x):
    return int((x+arb(1)/2).floor().unique_fmpz())


def simultaneous_approximants(alpha, beta, bound):
    proposals = set()
    independent = None
    for extra in range(0, 94, 6):
        k = ((1 << (10+extra))*(bound+3000))**3
        aa, bb = rounded(k*alpha), rounded(k*beta)
        reduced = fmpz_mat([[1, aa, bb], [0, k, 0], [0, 0, k]]).lll()
        rows = []
        for denominator, _, _ in reduced.tolist():
            d = abs(int(denominator))
            if not d:
                continue
            u, v = rounded(d*alpha), rounded(d*beta)
            if (128*(bound+3000)*abs(d*alpha-u) < 1
                    and 128*(bound+3000)*abs(d*beta-v) < 1):
                proposals.add(d)
                rows.append((u, v, d))
        if independent is None and len(rows) == 3 and fmpz_mat(rows).det() != 0:
            independent = max(row[2] for row in rows)
        if extra >= 36 and independent is not None:
            break
    assert independent is not None
    return sorted(proposals), independent


def divide_out(n, p):
    e = 0
    while n % p == 0:
        n //= p
        e += 1
    return n, e


def check_pass(p, q, filters, b, depth, following_bound):
    bound = 10**b
    ctx.prec = 3*bound.bit_length()+510
    log7 = arb(7).log()
    alpha, beta = arb(p).log()/log7, arb(q).log()/log7
    choices, homogeneous = simultaneous_approximants(alpha, beta, bound)
    powers = [1]
    while powers[-1] < 16384*max(choices):
        powers.append(7*powers[-1])
    assert 16384*max(choices) < 1 << 3000
    counts = [0]*(depth+1)
    largest_gap = 0
    log6 = arb(6).log()

    def visit(length, coefficient, gap):
        nonlocal largest_gap
        if any(not flags[length][coefficient % modulus] for modulus, flags in filters):
            return
        counts[length] += 1
        # j < 6 gives R = j mod 6. Our p,q are both primes above seven
        # in every nontrivial pair, so only the last level can be homogeneous.
        is_homogeneous = False
        if coefficient % 6 == 0:
            reduced, up = divide_out(coefficient//6, p)
            reduced, uq = divide_out(reduced, q)
            reduced, _ = divide_out(reduced, 7)
            is_homogeneous = reduced == 1
            if is_homogeneous:
                assert max(up, uq) <= 3000
        if is_homogeneous:
            denom = homogeneous
        else:
            shift = (arb(coefficient).log()-log6)/log7
            for denom in choices:
                value = denom*shift
                nearest = (value+arb(1)/2).floor().unique_fmpz()
                if nearest is not None and abs(value-int(nearest)) > arb(1)/32:
                    break
            else:
                raise ArithmeticError(("unseparated", p, q, b, length, coefficient))
        stop = bisect_left(powers, 16384*denom)
        if length == depth:
            largest_gap = max(largest_gap, stop)
        else:
            for next_gap in range(gap, stop):
                visit(length+1, coefficient*powers[next_gap-gap]+1, next_gap)

    visit(1, 1, 0)
    if following_bound is not None:
        c = 10**16*(2+3*b)
        assert (largest_gap+2)*c**(6-depth) < 10**following_bound
    return sum(counts), largest_gap


def small_weights(p, q):
    # E = bit length of 6t.
    limit = (1 << 3000)//6
    pa = p
    counts = [0, 0, 0]
    while pa*q <= limit:
        t = pa*q
        while t <= limit:
            counts[0] += 1
            if t % 504 == 83:
                counts[1] += 1
                if (6*t) & (t+1) == 0:
                    counts[2] += 1
                    n = 6*(t+1)
                    digit_sum = 0
                    while n and digit_sum <= 6:
                        n, digit = divmod(n, 7)
                        digit_sum += digit
                    assert digit_sum > 6
            t *= q
        pa *= p
    return counts


def base_pairs():
    # A repeated base or one equal to 7 leaves at most two predecessor primes.
    # A factor 3 is incompatible with t = 11 mod 72.
    primes = [p for p in primes_through(47) if p not in (2, 3, 7)]
    assert all(pow(p, 6, 504) == 1 for p in primes)
    return [(p, q) for p, q in combinations(primes, 2)
            if any(pow(p, a, 504)*pow(q, b, 504) % 504 == 83
                   for a in range(6) for b in range(6))]


if __name__ == "__main__":
    pairs = base_pairs()
    assert pairs == [(13, 23), (17, 43), (29, 31), (37, 47), (41, 43)]
    assert 2*(10**16*392)**6 < 10**130
    for p, q in pairs:
        filters = possible_prefixes(p, q)
        result = [check_pass(p, q, filters, *parameters) for parameters in PASSES]
        print((p, q), "prefix counts and gap cutoffs", result,
              "finite checks", small_weights(p, q), flush=True)
