"""Verify the two residual progressions at t=2661 by extracting quotient bits."""

from argparse import ArgumentParser
from math import isqrt


def check(limit):
    t, ratio = 2661, 42
    scale, stride = t*ratio, 10*t*ratio
    starts = (111761,357461)
    # Both possible CRT classes from x=2661 and x=3.
    for x,start in zip((2661,3),starts):
        assert start % (42*x) == 42*x-1
        assert (start-t) % (20*(t//x)) == 0
    popcount = getattr(int, "bit_count", lambda number: bin(number).count("1"))
    bits = 40
    mask = (1 << bits)-1
    checked = low_pass = 0
    for parity in range(2):
        top = limit-(limit-parity)%2
        if top < parity:
            continue
        value = scale << top
        prefix = isqrt(value//stride**2)
        assert (prefix*stride)**2 <= value < ((prefix+1)*stride)**2
        packed = prefix.to_bytes((prefix.bit_length()+7)//8,'little')
        for power in range(parity,top+1,2):
            shift = (top-power)//2
            window = int.from_bytes(packed[shift//8:shift//8+6],'little')
            low = (window >> (shift%8)) & mask
            for start in starts:
                for offset in (0,1):
                    e = stride*(low+offset)+start
                    a, rem = divmod(e*(e+1),scale)
                    assert rem == 0
                    u,v,w = a & mask,(a-e) & mask,(40*a) & mask
                    checked += 1
                    if popcount(u)+popcount(v)+popcount(w) != popcount(u+v+w):
                        continue
                    low_pass += 1
                    e = stride*((prefix >> shift)+offset)+start
                    a = e*(e+1)//scale
                    P = 1 << power
                    if a-e < P <= a < 2*P and (a | (a-e) | (40*a)) == 42*a-e:
                        raise ArithmeticError(('requires further tests',power,e))
    print('Broader quotient windows:',checked,'points,',low_pass,'passing the low binary test')


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('limit',type=int,nargs='?',default=1000000)
    check(parser.parse_args().limit)
