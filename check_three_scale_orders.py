"""Check digit exponents for exceptional orders at the pair cofactors 5, 7 and 9."""

from flint import fmpz

from carry_blocks import CarryGraph
from check_paired_power_orders import digit_matrix


def check_graph(base, slopes, offsets, numerator):
    length = 4096
    count = CarryGraph(base, slopes, offsets).block_bound(length)
    if len(slopes) == 2:
        matrix = digit_matrix(base, 1, *slopes)
        assert count == sum((matrix**length).entries())
    assert fmpz(count)**10000 <= fmpz(base)**(length*numerator)


def check():
    # Upper exponents in units of 1/10000, from unrestricted internal words.
    rows = (
        (5, 3, (2, 3, 5), (5979, 6788, 5033)),
        (7, 3, (2, 3, 7), (5522, 6472, 6109)),
        (7, 5, (2, 5, 7), (6311, 7281, 5807)),
        (9, 5, (2, 3, 5), (5737, 5453, 6982)),
        (9, 7, (2, 3, 7), (5647, 4840, 7573)),
        (9, 8, (2, 3), (6951, 68)),
        (9, 4, (2, 3), (5737, 5453)),
    )
    for d, k, bases, numerators in rows:
        for base, numerator in zip(bases, numerators):
            check_graph(base, (k, d-k), (0, 1), numerator)
        if len(bases) == 3:
            tau = max(sum(numerators)-10000, *numerators)
            threshold = {5: 9270, 7: 9800, 9: 9400}[d]
            assert 20000 + tau < 3*threshold
        print(d, k, bases, numerators)

    check_graph(2, (1, 2, 4), (0, 0, 0), 5519)
    assert 20000 + max(5519+4840+7573-10000, 7573) < 3*9400
    # The upper function (r+3)/(4 r^s) is increasing when r > 3s/(1-s).
    for limit, numerator, order_bound in (
        (2000, 8200, 9800), (100000, 8800, 9997), (110000, 8806, 9999)
    ):
        for r in (5, 7, 11, 13, 17, 19, 23, limit):
            choices = 1 + (r-1)//4 if r != limit else None
            if choices is None:
                assert fmpz(r+3)**10000 <= fmpz(4)**10000 * fmpz(r)**numerator
            else:
                assert fmpz(choices)**10000 <= fmpz(r)**numerator
        assert 23*(10000-numerator) >= 3*numerator
        tau = max(5737+5453+numerator-10000, numerator)
        assert 20000 + tau < 3*order_bound
        print('Prime divisor through', limit, 'digit exponent', numerator,
              'order exponent', order_bound)


if __name__ == '__main__':
    check()
