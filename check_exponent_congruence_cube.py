"""Prove the periodic extensions of the cube using the forced part and integer floors."""
import sys

BINARY_BITS = 64
TERNARY_DIGITS = 32
QUINARY_DIGITS = 24


def fits_floor(part, total, base, modulus):
    # Both arguments can be reduced independently to representatives.
    place = base
    while place <= modulus:
        if total//place != part//place + (total-part)//place:
            return False
        place *= base
    return True


def check_cube(limit):
    binary_mask = (1 << BINARY_BITS)-1
    ternary = 3**TERNARY_DIGITS
    quinary = 5**QUINARY_DIGITS
    powers_three = [3**c for c in range(1, limit+1)]
    powers_two = [2**b for b in range(limit+1)]
    binary_count = 0
    last_cases, final = [], []
    D = 1
    for a in range(1, limit+1):
        D *= 5
        total3 = [D*Q % ternary for Q in powers_two]
        inverse = pow(2, -limit, D)
        third_inverse = pow(3, -1, D)
        for c, R in enumerate(powers_three, 1):
            inverse = inverse*third_inverse % D
            total = D*R
            forbidden = (~total) & binary_mask
            u = inverse*R
            for b in range(limit, 0, -1):
                if u & forbidden == 0:
                    x = u//R
                    binary_count += 1
                    Q = powers_two[b]
                    if fits_floor(x*Q % ternary, total3[b], 3, ternary):
                        W = Q*R
                        quotient, remainder = divmod(x*W-1, D)
                        assert remainder == 0
                        if fits_floor(quotient % quinary, W % quinary, 5, quinary):
                            final.append((a, b, c))
                        else:
                            depth = next(j for j in range(1, QUINARY_DIGITS+1)
                                         if quotient % 5**j > W % 5**j)
                            last_cases.append((a, b, c, depth))
                u += u
                if u >= total:
                    u -= total
    return binary_count, sorted(last_cases), final


if __name__ == '__main__':
    limit = int(sys.argv[1]) if len(sys.argv) > 1 else 1024
    assert 1 <= limit <= 1024
    count, cases, final = check_cube(limit)
    print(limit**3, 'cells;', count, 'passing binary;', len(cases)+len(final), 'passing both')
    print('Base-five cases (a, b, c, depth):', cases)
    assert not final, final
    if limit == 1024:
        assert count == 3182 and len(cases) == 7
        assert max(a+j for a, b, c, j in cases) == 9
