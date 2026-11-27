"""Decide bounded prime valuation near any positive rational trinomial direction."""
import argparse
from math import gcd

def bounded_path(slopes, p, budget, *, witness_caps=None):
    """Decide bounded valuation. Optional tighter caps only search for witnesses."""
    total = sum(slopes)
    lower = [a-total for a in slopes]
    # At zero budget every unread tail is carry-free, so its total is at
    # most one. This gives the sharper upper discrepancy bound s_i.
    upper = list(slopes) if budget == 0 else [2*a for a in slopes]
    if witness_caps is not None:
        upper = [min(a,b) for a,b in zip(upper,witness_caps)]

    def children(state, first=False):
        z1,z2,out,remaining = state
        if remaining < out:
            return
        z = (z1,z2,-z1-z2)
        for digit in range(int(first),p):
            for incoming in range(min(2,remaining-out)+1):
                column = p*out+digit-incoming
                if not 0 <= column <= 3*(p-1):
                    continue
                limits = []
                for i in (0,1):
                    base = p*z[i] - slopes[i]*column
                    limits.append(range(max(0,-((base-lower[i])//total)),
                                        min(p-1,(upper[i]-base)//total)+1))
                for a in limits[0]:
                    for b in limits[1]:
                        c = column-a-b
                        if not 0 <= c < p:
                            continue
                        delta = [p*old+total*x-s*column
                                 for old,x,s in zip(z,(a,b,c),slopes)]
                        if lower[2] <= delta[2] <= upper[2]:
                            yield (delta[0],delta[1],incoming,remaining-out),(a,b,c)

    # A separate root requires a nonzero first digit of the total even when
    # that digit comes entirely from the incoming carry.
    origin = (0,0,0,budget)
    states = [None]
    edges = [iter(children(origin,True))]
    digit_path = []
    active = {}
    done = set()
    while states:
        edge = next(edges[-1],None)
        if edge is None:
            old = states.pop()
            edges.pop()
            if old is not None:
                del active[old]
                done.add(old)
                digit_path.pop()
            continue
        target, digits = edge
        if target in done:
            continue
        if target in active:
            path = digit_path+[digits]
            cut = active[target]
            return path[:cut], path[cut:]
        active[target] = len(states)
        states.append(target)
        edges.append(iter(children(target)))
        digit_path.append(digits)
    return None


def floor_value(parts, p):
    value = 0
    m = p
    total = sum(parts)
    while m <= total:
        value += total//m - sum(a//m for a in parts)
        m *= p
    return value


def check_period(slopes,p,budget,witness):
    prefix, cycle = witness
    def row(digits):
        parts = [0,0,0]
        for word in digits:
            parts = [p*v+d for v,d in zip(parts,word)]
        return parts
    initial, block = row(prefix), row(cycle)
    modulus = p**len(cycle)-1
    vector = [modulus*a+b for a,b in zip(initial,block)]
    assert vector[0] > 0 and all(x*slopes[0] == vector[0]*s for x,s in zip(vector,slopes))
    repeat = 1
    for _ in range(3):
        actual = [a+repeat*b for a,b in zip(initial,vector)]
        assert all(actual) and floor_value(actual,p) <= budget
        repeat = (modulus+1)*repeat+1


def forced_divisor(slopes):
    from check_equal_part_cones import is_prime
    total = sum(slopes)
    result = 1
    for p in range(2,total+1):
        if not is_prime(p):
            continue
        upper = floor_value(slopes,p)
        n = total
        while n%p == 0:
            upper -= 1
            n //= p
        assert upper >= 0  # The primitive sum divides the coefficient.
        for budget in range(upper):
            witness = bounded_path(slopes,p,budget)
            if witness:
                check_period(slopes,p,budget,witness)
                result *= p**budget
                break
        else:
            result *= p**upper
    return result


def cross_checks():
    from check_equal_part_cones import difference_of_powers, graph_witness, is_prime
    tests = 0
    for m in range(1,251):
        for p in range(2,m+2):
            if (m+1)%p or not is_prime(p):
                continue
            n = m+1
            v = 0
            while n%p == 0:
                n //= p
                v += 1
            forced = p == 2 or difference_of_powers(m+1,p)
            assert forced == (graph_witness(m,p) is None)
            for budget in {0,v-1,v}:
                witness = bounded_path((1,1,m),p,budget)
                assert (witness is None) == (forced and budget<v), (m,p,budget)
                if witness:
                    check_period((1,1,m),p,budget,witness)
                tests += 1
    for slopes in ((1,2,4),(1,3,7),(5,13,19),(2,3,5),(4,5,9)):
        for p in (2,3,5,7):
            for budget in range(floor_value(slopes,p)+1):
                witness = bounded_path(slopes,p,budget)
                if witness:
                    check_period(slopes,p,budget,witness)
                    break
            else:
                raise AssertionError((slopes,p))
    for slopes,expected in (((1,2,4),3),((1,3,7),40),((5,13,19),14)):
        assert forced_divisor(slopes) == expected
    print('Independent proportion-graph comparisons:',tests)


def cofactor_five_controls():
    from math import isqrt
    for q,z in ((1467779,2918400),(2823523,5636100)):
        slopes = (3*q,z,2*q-z)
        lengths = []
        for p in (2,3,5):
            a,b,c = slopes
            witness = bounded_path(slopes,p,0,witness_caps=(a-1,b-1,c))
            assert witness
            check_period(slopes,p,0,witness)
            lengths.append(tuple(map(len,witness)))
        assert all(q%r for r in range(2,isqrt(q)+1))
        assert floor_value((3*q,z-1,2*q-z),q) == 0
        print('Cofactor-five prime and prefix/period lengths in bases 2,3,5:',q,lengths)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('parameters',type=int,nargs='*',metavar='integer',
                        help='a b c for the divisor, optionally followed by prime and maximum valuation')
    args = parser.parse_args()
    if args.parameters:
        from check_equal_part_cones import is_prime
        if len(args.parameters) not in (3,5):
            parser.error('provide a b c, optionally followed by prime and maximum valuation')
        a,b,c = args.parameters[:3]
        if min(a,b,c)<1:
            parser.error('positive slopes are required')
        g = gcd(a,b,c)
        slopes = (a//g,b//g,c//g)
        if len(args.parameters) == 3:
            print('Largest divisor forced in a neighborhood:',forced_divisor(slopes))
            raise SystemExit(0)
        p,budget = args.parameters[3:]
        if budget<0 or not is_prime(p):
            parser.error('a prime base and nonnegative valuation are required')
        result = bounded_path(slopes,p,budget)
        if result is None:
            print(f'A neighborhood forces divisibility by {p}^{budget+1}.')
        else:
            check_period(slopes,p,budget,result)
            print('Prefix and repeating block of part digits:',result)
    else:
        cross_checks()
        cofactor_five_controls()
