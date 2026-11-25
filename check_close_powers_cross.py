"""Cross-check square-case assignments by weight and complementary rows by part."""

from math import gcd, isqrt, lcm

from check_close_powers import INITIAL_LIMIT, square
from check_six_bit_rays import normalized_rows
from check_six_digit_gaps import complement_by_parts
from prime_certificates import primes_through


def divisors_by_trial(n):
    small = [k for k in range(1, isqrt(n)+1) if n % k == 0]
    return small + [n//k for k in reversed(small) if k*k != n]


def count_at_scale(P, Q, v, layers):
    d = Q-P
    predecessor = P*Q-1
    lower = isqrt((27*predecessor**2-1)//(P*v)**3)+1
    upper = d*d//(4*v*v)
    if upper < lower:
        return 0
    available = gcd(predecessor, layers[d//v])
    if available < lower:
        return 0
    count = 0
    for s in divisors_by_trial(available):
        if not lower <= s <= upper:
            continue
        parts = [x for x in range(v, d, v) if x*(x-Q) % s == 0]
        for i, x in enumerate(parts):
            for y in parts[i:]:
                if x+y > d:
                    break
                z = Q-x-y
                if lcm(gcd(x,y), gcd(x,z), gcd(y,z)) % s == 0:
                    count += 1
    return count


def weight_first(max_gap):
    cap = max_gap**4//432
    layers = [1]
    for j in range(1, max_gap+1):
        layers.append(lcm(j, layers[-1]))
    count = 0
    for p in primes_through(cap):
        P = p
        while P <= cap:
            if P*(P+max_gap) > INITIAL_LIMIT:
                # Combining the complementary size bound with 4*v**2*s <= d**2
                # gives 432*P*v < d**4 since P*(P+d)-1 > P**2.
                for v in range(1, min(max_gap//2, cap//P)+1):
                    if v % p == 0:
                        continue
                    start = max(2*v, isqrt(isqrt(432*P*v))+1)
                    for d in range(start+(-P-start) % v, max_gap+1, v):
                        Q = P+d
                        if d % p and P*Q > INITIAL_LIMIT and gcd(Q//v, v) == 1:
                            count += count_at_scale(P, Q, v, layers)
            P *= p
    return count


def compare_complements(total, powers):
    by_part = complement_by_parts(total, powers)
    by_position = {tuple(sorted(row)) for row in normalized_rows(total, powers)}
    if by_part != by_position:
        raise ArithmeticError(("complement mismatch", total, powers))
    return by_part


def positive_controls():
    for row in [(20, 36, 45), (200, 216, 225)]:
        total = sum(row)
        powers = (4, 5, 9)
        direct = set()
        for a in range(1, total//3+1):
            for b in range(a, (total-a)//2+1):
                parts = a, b, total-a-b
                if all(sum(x % w == 0 for x in parts) == 2 for w in powers):
                    direct.add(parts)
        if row not in direct or compare_complements(total, powers) != direct:
            raise ArithmeticError(("positive control", row))


if __name__ == "__main__":
    positive_controls()
    reverse = weight_first(160)
    forward = square(160)
    if reverse != forward["full prime-power assignments"]:
        raise ArithmeticError(("weight-first mismatch", reverse, forward))
    print(f"Weight-first enumeration agrees at {reverse} assignments.")
    square(300, compare_complements)
