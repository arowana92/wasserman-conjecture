"""Exclude huge binary exponent intervals by a rational grid and ternary digit words."""

from argparse import ArgumentParser
from fractions import Fraction

from check_seven_quadratic_bases_cross import S, ceil_div, log_interval
from check_seven_quadratic_suffixes import RESIDUAL
from check_seven_ternary_words import graph_data, leading_graph

F = 3**220
CHECK_TRITS = 125
LIMITS = {4587:[Fraction(256,417)],7383:[Fraction(19597,35328),Fraction(486233,876544),
          Fraction(836139,1507328),Fraction(46678369,84148224)],9453:[Fraction(256,411)]}


def natural(n,d):
    """Use the rational logarithm enclosure on ratios in [1,3]."""
    if n>2*d:
        lo,hi = log_interval(2,1)
        d *= 2
    else:
        lo = hi = 0
    k = 256*n//d
    x,y = log_interval(k,256)
    u,v = log_interval(256*n,k*d)
    return lo+x+u,hi+y+v


LOG3 = natural(3,1)


def log3(n,d):
    exponent = 0
    while n >= 3*d:
        d *= 3
        exponent += 1
    while n < d:
        n *= 3
        exponent -= 1
    x,y = natural(n,d)
    return exponent*S+x*S//LOG3[1], exponent*S+ceil_div(y*S,LOG3[0])


def approximation(bound):
    lo,hi = log3(2,1)
    x,y = lo,S
    prev_p,p,prev_q,q = 0,1,1,0
    while y:
        a,r = divmod(x,y)
        prev_p,p,prev_q,q = p,a*p+prev_p,q,a*q+prev_q
        if q>bound and q*max(abs(q*lo-p*S),abs(q*hi-p*S)) <= S:
            return p,q
        x,y = y,r
    raise ValueError('insufficient precision for the exponent grid')


def words_at(first,transitions,size):
    current = first
    for _ in range(1,size):
        current = [(3*word+d,dst) for word,state in current for d,dst in transitions[state]]
    return current


def integer_grid_candidates(words,depth,phase,grid,bound):
    beta_lo,beta_hi = phase
    a,q = grid
    inv = pow(a,-1,q)
    base = 3**(depth-1)
    candidates = set()
    for word,_ in words:
        lo = log3(word,base)[0]
        hi = log3(word+1,base)[1]
        start = (q*(lo-beta_hi))//S-2
        stop = ceil_div(q*(hi-beta_lo),S)+2
        for i in range(start,stop+1):
            m = inv*i%q
            if 65536 <= m <= bound:
                candidates.add(m)
    return candidates


def normalize(lo,hi):
    while lo>=3*F:
        lo //= 3
        hi = ceil_div(hi,3)
    return lo,hi


def multiply(first,second):
    return normalize(first[0]*second[0]//F, ceil_div(first[1]*second[1],F))


def powers_of_two(bits):
    word = (2*F,2*F)
    table = []
    for _ in range(bits):
        table.append(word)
        word = multiply(word,word)
    return table


def mantissa(m,scale,table):
    x = scale.numerator*F
    word = normalize(x//scale.denominator,ceil_div(x,scale.denominator))
    while m:
        bit = (m & -m).bit_length()-1
        word = multiply(word,table[bit])
        m &= m-1
    return word


def mantissa_intervals(lo,hi):
    # The projected true scale differs by less than a single unit of this precision.
    lo -= 1
    hi += 1
    result = [(max(F,lo),min(3*F,hi))]
    if lo < F:
        result.append((3*lo,3*F))
    if hi >= 3*F:
        result.append((F,ceil_div(hi,3)))
    return [(x,y) for x,y in result if x<=y]


def accepted(word,depth,first,arcs):
    digits = []
    for _ in range(depth):
        word,d = divmod(word,3)
        digits.append(d)
    if word or digits[-1] not in first:
        return False
    state = first[digits.pop()]
    for d in reversed(digits):
        if d not in arcs[state]:
            return False
        state = arcs[state][d]
    return True


def leading_possible(value,first,arcs):
    place = 3**(CHECK_TRITS-1)
    for lo,hi in mantissa_intervals(*value):
        a,b = max(place,lo*place//F),min(3*place-1,hi*place//F)
        assert b-a<=2
        if any(accepted(w,CHECK_TRITS,first,arcs) for w in range(a,b+1)):
            return True
    return False


def run(decimal_digits):
    bound = 2*10**decimal_digits
    grid = approximation(bound)
    depth = 1
    while 3**(depth-1)<grid[1]:
        depth += 1
    table = powers_of_two(bound.bit_length())
    assert grid[1]<10**34 and CHECK_TRITS < 40000
    print(f'm <= {bound}, grid {grid[0]}/{grid[1]}, first cut {depth} ternary digits.',flush=True)
    for D,_,_,z,w in RESIDUAL:
        back,terminal,_ = graph_data(D,z,w)
        starts,transitions = leading_graph(back,terminal)
        words = words_at(starts,transitions,depth)
        first,arcs = dict(starts),list(map(dict,transitions))
        print(f'D={D}: {len(words)} initial words.',flush=True)
        for f in LIMITS[D]:
            scale = 16*f
            phase = log3(scale.numerator,scale.denominator)
            indices = integer_grid_candidates(words,depth,phase,grid,bound)
            left = [m for m in indices if leading_possible(mantissa(m,scale,table),first,arcs)]
            print(f'  limit {f}: {len(indices)} grid candidates, {len(left)} at {CHECK_TRITS} trits.',flush=True)
            assert not left, left[:5]


if __name__ == '__main__':
    cli = ArgumentParser(description=__doc__)
    cli.add_argument('digits',type=int,nargs='?',default=30)
    args = cli.parse_args()
    if not 6 <= args.digits <= 30:
        cli.error('use an exponent limit from 10^6 through 10^30')
    run(args.digits)
