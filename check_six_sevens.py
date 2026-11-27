"""Prove the cutoff for six positive powers of seven near a power of two."""
from bisect import bisect_left, bisect_right
from math import gcd, isqrt

from flint import arb, ctx


# (initial upper bound for E, number of prefixes to separate, new bound for E)
ROUNDS = ((10**180, 2, 10**74), (10**74, 3, 10**55), (10**55, 4, 10**38),
          (10**38, 5, 10**21), (10**21, 6, 3000), (3000, 6, 200))


def approximants(alpha, bound):
    # A dyadic midpoint only suggests fractions. Each accepted fraction is
    # checked for coprimality, size and accuracy against the original ball.
    denominator = 1 << (ctx.prec-30)
    numerator = int((alpha*denominator).floor().unique_fmpz())
    h, q, old_h, old_q = 1, 0, 0, 1
    choices = []
    while denominator:
        partial, remainder = divmod(numerator, denominator)
        h, old_h = partial*h+old_h, h
        q, old_q = partial*q+old_q, q
        if q > bound*10**27:
            break
        if q > 1000*(bound+3000) and 64*(bound+3000)*abs(q*alpha-h) < 1:
            assert gcd(q, h) == 1
            choices.append(q)
        numerator, denominator = denominator, remainder
    assert choices
    return choices


def check_round(bound, depth, target):
    ctx.prec = 2*bound.bit_length()+240
    log7 = arb(7).log()
    alpha = arb(2).log()/log7
    choices = approximants(alpha, bound)
    powers = [1]
    while powers[-1] < 12800*choices[-1]:
        powers.append(7*powers[-1])
    counts = [0]*(depth+1)
    gaps = [0]*(depth+1)
    largest_q = 0
    # Throughout the first five passes E > 3000 suffices. The last assumes E > 200.
    threshold = 200 if bound == 3000 else 3000

    def visit(length, coefficient, gap):
        nonlocal largest_q
        counts[length] += 1
        reduced = coefficient//(coefficient & -coefficient)
        while reduced % 7 == 0:
            reduced //= 7
        if reduced == 1:
            two_exponent = (coefficient & -coefficient).bit_length()-1
            assert two_exponent <= threshold
            q = choices[0]
        else:
            shift = arb(coefficient).log()/log7
            for q in choices:
                value = q*shift
                nearest = (value+arb(1)/2).floor().unique_fmpz()
                if nearest is not None and abs(value-int(nearest)) > arb(1)/32:
                    break
            else:
                raise ArithmeticError(("no separation certificate", coefficient))
        largest_q = max(largest_q, q)
        assert 12800*q < 1 << (threshold//2)
        stop = bisect_left(powers, 12800*q)
        gaps[length] = max(gaps[length], stop)
        if length < depth:
            for next_gap in range(gap, stop):
                new_coefficient = coefficient*powers[next_gap-gap]+1
                visit(length+1, new_coefficient, next_gap)

    visit(1, 1, 0)
    if depth < 6:
        # C = 10^14(1+log(E+1)). Upper-bound log(10) by 3.
        c = 10**14*(2+3*(len(str(bound))-1))
        assert (gaps[depth]+2)*c**(6-depth) < target
    else:
        assert threshold == target
    print(f"E <= {bound:.1e}, depth {depth}: {sum(counts):,} prefixes, "
          f"gap cutoffs {gaps[1:]}, denominator at most {largest_q.bit_length()} bits.")
    return sum(counts)


def exact_hits(exponent):
    target = 2**exponent
    distance = isqrt(16*target)
    lower, upper = target-distance, target+distance
    powers = []
    power = 7
    while power <= upper:
        powers.append(power)
        power *= 7

    def search(total, slots, stop):
        if lower <= total <= upper and total:
            yield total
        if slots and stop:
            smallest = (lower-total+slots-1)//slots
            for j in range(bisect_left(powers, smallest),
                           bisect_right(powers, upper-total, 0, stop)):
                yield from search(total+powers[j], slots-1, j+1)
    return list(search(0, 6, len(powers)))


if __name__ == "__main__":
    assert 2*(10**14*542)**6 < 10**180
    checked = sum(check_round(*row) for row in ROUNDS)
    assert all(not exact_hits(e) for e in range(18, 201))
    assert exact_hits(17) == [129654]
    print(f"{checked:,} certified prefixes. The last hit is 129654 at exponent 17.")
