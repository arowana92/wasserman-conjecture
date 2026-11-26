"""Exact binary-power windows for the first remaining resonant gap multiplier."""

from argparse import ArgumentParser
from collections import Counter
from math import gcd, isqrt

from check_resonant_row_gaps_cross import audit, digits_pass


def binary_windows(params, limit, direct):
    x,y,b,c,m,h,n,j = params
    t = x*y
    scale = t*(c//n)
    assert c % n == 0 and j % (2*n) == 0
    # e+1 = xc*T and e-t = yb*S.
    # Equivalently e is -1 modulo xc and t modulo yb.
    left,right = x*c,y*b
    common = gcd(left,right)
    step = (left//common)*right
    multiplier = ((t+1)//common)*pow(left//common,-1,right//common)
    origin = (left*multiplier-1) % step
    assert (origin+1) % left == (origin-t) % right == 0
    layer = 1 << 40
    modulus = (step//gcd(step,scale*layer))*(scale*layer)
    counts = Counter()

    for parity in range(2):
        top = limit-(limit-parity)%2
        if top < parity:
            continue
        high_root = isqrt(scale << top)
        assert high_root**2 <= scale << top < (high_root+1)**2
        packed = high_root.to_bytes((high_root.bit_length()+7)//8,'little')
        residue = isqrt(scale << parity) % modulus
        for exponent in range(parity,top+1,2):
            shift = (top-exponent)//2
            if direct:
                residue = isqrt(scale << exponent)
            elif exponent != parity:
                bit = (packed[shift//8] >> (shift%8)) & 1
                residue = (2*residue+bit) % modulus
            for offset in range((origin-residue)%step,scale+1,step):
                counts['progression points in windows'] += 1
                reduced = (residue+offset) % (scale*layer)
                leading,rem = divmod(reduced*(reduced+1),scale)
                assert rem == 0
                row = leading,leading-reduced,(c//n-2)*leading
                if direct:
                    passes = digits_pass(row,2,layer)
                else:
                    a0,a1,a2 = (v % layer for v in row)
                    passes = not ((a0 & a1) | (a1 & a2) | (a2 & a0))
                if not passes:
                    continue
                counts['passing 40 binary digits'] += 1
                e = (high_root >> shift)+offset
                large = e*(e+1)//scale
                P = 1 << exponent
                if e <= t or not large-e < P <= large < 2*P:
                    continue
                actual = large,large-e,(c//n-2)*large
                if all(k > 0 for k in actual) and digits_pass(actual,2,1 << sum(actual).bit_length()):
                    raise ArithmeticError(('requires remaining carry tests',params,exponent,e))
    print('x, y',x,y,'through exponent',limit,dict(counts))
    return counts


def main(limit, direct):
    if limit < 0:
        raise ValueError('use a nonnegative exponent bound')
    def check(params,bins):
        binary_windows(params,limit,direct)
    audit(2661,first_t=2661,check_remaining=check)


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('limit', type=int, nargs='?', default=1000000)
    parser.add_argument('--direct',action='store_true',help='compute a fresh integer square root at every exponent')
    args = parser.parse_args()
    main(args.limit,args.direct)
