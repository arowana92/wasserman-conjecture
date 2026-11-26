"""Independent factorial-floor and rational-reconstruction check at 986."""
from argparse import ArgumentParser
from itertools import combinations
from math import gcd


def layer_is_free(slopes, x, modulus):
    a, b, c = slopes
    return a*x % modulus + b*x % modulus + (c*x+1) % modulus < modulus


def prefixes(slopes, k):
    values = [1]
    modulus = 7
    for _ in range(1,k):
        values = [z for x in values for z in range(x,x+7*modulus,modulus)
                  if layer_is_free(slopes,z,7*modulus)]
        modulus *= 7
    return values


def log8(x, k):
    exponent = 0
    step = 1
    for _ in range(1,k):
        modulus = step*49
        digit = (x-pow(8,exponent,modulus))//(step*7) % 7
        exponent += digit*step
        step *= 7
    assert pow(8,exponent,7**k) == x
    return exponent


def reconstruct(residue, modulus, bound):
    """Recover a possible positive reduced numerator and denominator <= bound."""
    upper, lower = modulus, residue
    old, coefficient = 0, 1
    while lower > bound:
        quotient, remainder = divmod(upper,lower)
        upper, lower = lower, remainder
        old, coefficient = coefficient, old-quotient*coefficient
    if lower and 0 < coefficient <= bound and gcd(lower,coefficient) == 1:
        return lower, coefficient
    return None


def lift(slopes, x, y, modulus, a, b):
    """Require x**b == y**a; compute the possible next y digit from the x digit."""
    nodes = [(x,y)]
    levels = 0
    inverse = pow(a,-1,7)
    while nodes:
        next_modulus = 7*modulus
        following = []
        for x,y in nodes:
            gap = (pow(x,b,next_modulus)-pow(y,a,next_modulus)) % next_modulus
            assert gap % modulus == 0
            diff = gap//modulus
            for i in range(7):
                j = (diff+b*i)*inverse % 7
                xx, yy = x+i*modulus, y+j*modulus
                if layer_is_free(slopes,xx,next_modulus) and layer_is_free(slopes,yy,next_modulus):
                    following.append((xx,yy))
        nodes = following
        levels += 1
        modulus = next_modulus
        assert levels <= 64, ("unresolved exponents",a,b)
    return levels


def verify(slopes, bound):
    k = 2
    while 68*7**(k-1) <= 2*bound**2:
        k += 1
    period = 7**(k-1)
    modulus = 68*period
    inverse = pow(period,-1,68)
    points = []
    for x in prefixes(slopes,k):
        z = log8(x,k)
        assert z%7 != 0
        points.append((x,z,pow(z,-1,period)))
    count = depth = 0
    for (x,z,_),(y,_,iz) in combinations(points,2):
        residue = z*iz % period
        residue += period*((1-residue)*inverse % 68)
        pair = reconstruct(residue,modulus,bound)
        if pair:
            a, b = pair
            assert a != b and (a-residue*b)%modulus == 0
            count += 1
            depth = max(depth,lift(slopes,x,y,7*period,a,b))
    print(slopes,"initial digits",k,"floor-admissible words",len(points),
          "reconstructed pairs with prefixes",count,"maximum extra digits",depth,flush=True)


def extendible_prefix(slopes, x, p, k):
    m = p**k
    if any(not layer_is_free(slopes,x,p**j) for j in range(1,k+1)):
        return False
    a,b,c = slopes
    start = (a*x//m, b*x//m, (c*x+1)//m)
    stack = [start]
    seen = {start}
    while stack:
        state = stack.pop()
        if state == (0,0,0):
            return True
        for digit in range(p):
            values = tuple(v+digit*t for v,t in zip(state,slopes))
            if sum(w%p for w in values) >= p:
                continue
            target = tuple(w//p for w in values)
            if target not in seen:
                seen.add(target)
                stack.append(target)
    return False


def check_initial_residues():
    for slopes, first, seven in (((119,203,664),1,{22,29}), ((203,714,69),1,{22})):
        for p, low, k, expected in ((7,first,2,seven),(17,16,2,{84})):
            assert {x for x in range(low,p**k,p)
                    if extendible_prefix(slopes,x,p,k)} == expected
    assert {x for x in range(1,7**5,7)
            if extendible_prefix((203,714,69),x,7,5)} == {5608}


def controls():
    for bound in range(2,15):
        modulus = 2*bound**2+1
        expected = {a*pow(b,-1,modulus)%modulus:(a,b)
                    for a in range(1,bound+1) for b in range(1,bound+1)
                    if gcd(a,b)==gcd(b,modulus)==1}
        for r in range(modulus):
            assert reconstruct(r,modulus,bound) == expected.get(r)
    for k in range(1,5):
        for e in range(7**(k-1)):
            assert log8(pow(8,e,7**k),k) == e
    print("Direct small reconstruction and discrete-log controls passed.",flush=True)


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--limit',type=int,default=10**6)
    args = parser.parse_args()
    controls()
    check_initial_residues()
    for row in ((119,203,664),(203,714,69)):
        verify(row,args.limit)
