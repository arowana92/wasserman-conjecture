"""Independent exponent enumeration by exact counts of residues in angular intervals."""

from argparse import ArgumentParser
from fractions import Fraction
from functools import lru_cache
from itertools import product

from check_seven_quadratic_bases import log_floor
from check_seven_quadratic_suffixes import RESIDUAL
from check_seven_ternary_rotation import LIMITS
from check_seven_ternary_words import no_carry, offsets

BITS = 365
S = 1 << BITS
LOG_TWO_THREE = log_floor(Fraction(3),BITS)


def ceiling(a,b):
    return -(-a//b)


def coarse_log(x):
    """Bracket S*log_3(x) for positive rational x >= 1."""
    lo = log_floor(x,BITS)
    return S*lo//(LOG_TWO_THREE+1), ceiling(S*(lo+1),LOG_TWO_THREE)


@lru_cache(maxsize=None)
def grid_log(k):
    return coarse_log(Fraction(k,256))


def inverse_ln3():
    # ln(3) = sum_{j>=0} 1/(4^j*(2j+1)).
    base = 1
    lo = hi = 0
    for j in range(BITS):
        d = base*(2*j+1)
        lo += S//d
        hi += ceiling(S,d)
        base *= 4
        if base>S*4:
            break
    hi += 1
    return S*S//hi, ceiling(S*S,lo)


INVLOG = inverse_ln3()


def endpoint(n,d):
    k = 256*n//d
    gl,gu = grid_log(k)
    a,b = 256*n-k*d,k*d
    if not a:
        return gl,gu
    # Alternating Taylor series at 1 with 33 terms; 0 <= a/b < 1/256.
    p = S*a//b
    q = ceiling(S*a,b)
    lo = hi = 0
    for j in range(1,34):
        if j%2:
            lo += p//j
            hi += ceiling(q,j)
        else:
            lo -= ceiling(q,j)
            hi -= p//j
        p = p*a//b
        q = ceiling(q*a,b)
    lo -= ceiling(q,34)
    assert lo >= 0
    return gl+lo*INVLOG[0]//S, gu+ceiling(hi*INVLOG[1],S)


def reverse_graph(D,z,w):
    slopes,states = offsets(D,z,w)
    universe = set(states)
    terminal = frozenset(s for s in states if no_carry(s))
    outputs = [d for d in product(range(3),repeat=3) if sum(d)<3]
    def prepend(group,k):
        old = set()
        for x,y,v in group:
            a,b,c = 3*x-k*slopes[0],3*y-k*slopes[1],3*v-k*slopes[2]
            for e,f,g in outputs:
                s = (a+e,b+f,c+g)
                if s in universe:
                    old.add(s)
        return frozenset(old)

    groups,idx,edges = [],{},[]
    def save(s):
        if s not in idx:
            idx[s]=len(groups)
            groups.append(s)
        return idx[s]
    roots = [(k,save(s)) for k in (1,2) if (s := prepend(terminal,k))]
    for group in groups:
        edges.append([(k,save(s)) for k in range(3) if (s := prepend(group,k))])
    return roots,edges


def floor_sum(n,m,a,b):
    """Sum floor((a*j+b)/m) for 0 <= j < n by exchanging the lattice coordinates."""
    h,b = divmod(b,m)
    result = h*n
    while True:
        h,a = divmod(a,m)
        result += h*n*(n-1)//2
        h,b = divmod(b,m)
        result += h*n
        y = a*n+b
        if y<m:
            return result
        n,b = divmod(y,m)
        m,a = a,m


def cross(bound):
    alpha,upper = coarse_log(Fraction(2))
    start, n = 65536,bound-65536+1
    assert n>0
    for D,_,_,z,w in RESIDUAL:
        roots,edges = reverse_graph(D,z,w)
        for f in LIMITS[D]:
            beta,top = coarse_log(16*f)
            shift = alpha*start+beta
            # Any true phase differs modulo one from (m*alpha+beta)/S by
            # at most 'radius'/S, including the true quadratic displacement.
            radius = bound*(upper-alpha)+(top-beta)+1
            nodes = roots
            place = 1
            for h in range(1,126):
                kept = []
                for word,state in nodes:
                    left = endpoint(word,place)[0]-radius
                    right = endpoint(word+1,place)[1]+radius+1
                    # Wide coarse intervals are harmless and need no counting.
                    if (right-left)*n>200*S:
                        kept.append((word,state))
                        continue
                    hits = floor_sum(n,S,alpha,shift-left)-floor_sum(n,S,alpha,shift-right)
                    if hits:
                        kept.append((word,state))
                if not kept:
                    print(f'D={D}, limit {f}: every index excluded at prefix length {h}.', flush=True)
                    break
                nodes = [(3*word+d,dst) for word,state in kept for d,dst in edges[state]]
                place *= 3
            else:
                raise AssertionError((D,f,len(kept)))


if __name__ == '__main__':
    cli = ArgumentParser(description=__doc__)
    cli.add_argument('digits',type=int,nargs='?',default=30)
    args = cli.parse_args()
    assert 6<=args.digits<=30
    for size in (1,7,20):
        for a,b in ((3,4), (20,-9101), (0,-4), (112,83)):
            assert floor_sum(size,17,a,b) == sum((a*i+b)//17 for i in range(size))
    cross(2*10**args.digits)
