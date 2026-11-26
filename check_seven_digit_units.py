"""Check the finite bounds and the digit allocations in the three-position theorem."""

from collections import Counter
from itertools import permutations, product
from math import gcd, lcm, prod

from check_dense_offsets import free_rows, has_size_obstruction

# For each row give the multiplicities of L,K,A, the two unshifted parts
# as linear forms, and the signed vectors at U and the two predecessor gcds.
# None means the corresponding odd predecessor factor is one.
G1, G2, G3, G4 = (-1,1,1), (-1,-1,2), (-1,-3,3), (-2,4,1)
CASES = (
    ((3,1,2), (0,1,1),(1,0,1), ((-1,-1,1),G3,G1)),
    ((3,1,2), (0,1,1),(2,0,1), ((-1,-2,2),G3,G2)),
    ((3,1,2), (1,0,1),(2,1,0), ((-1,2,1),G1,G4)),
    ((3,1,2), (1,0,2),(1,1,0), ((-2,2,1),G4,G1)),
    ((3,1,2), (1,1,0),(2,0,1), ((-1,1,2),G1,G2)),
    ((4,1,1), (1,1,0),(1,0,1), ((1,-1,-1),(1,-1,-3),(1,-3,-1))),
    ((4,1,1), (1,1,0),(2,0,1), ((1,-1,-2),(1,-1,-3),(1,-2,-2))),
    ((2,2,2), (1,1,0),(0,1,2), ((2,-2,1),None,(1,-2,1))),
)


def dot(form, values):
    return sum(a*b for a,b in zip(form,values))


def unshifted_pairs(mult):
    kinds = [kind for kind, amount in enumerate(mult) for _ in range(amount)]
    found = set()
    for destinations in product(range(3), repeat=6):
        slots = [[0, 0, 0] for _ in range(3)]
        for kind, slot in zip(kinds, destinations):
            slots[slot][kind] += 1
        if not any(slots[2]) or any(sum(v > 0 for v in slot) < 2 for slot in slots[:2]):
            continue
        found.add(tuple(sorted(tuple(slot) for slot in slots[:2])))
    return found


def template_controls():
    classified = {mult: set() for mult in ((3,1,2),(4,1,1),(2,2,2))}
    for mult,x,y,points in CASES:
        classified[mult].add((x,y))
        for pair, point in zip(((x,y),(mult,x),(mult,y)),points):
            if point is not None:
                assert all(dot(form,point) == 0 for form in pair)
    classified[(3,1,2)].update((
        ((1,0,1),(1,0,1)), ((1,0,1),(1,1,0)),
        ((1,0,1),(1,1,1)), ((1,0,1),(2,0,1))))
    classified[(2,2,2)].update((
        ((1,1,0),(1,1,0)), ((1,1,0),(1,0,1)), ((1,1,0),(1,1,1))))
    for mult, pairs in classified.items():
        possibilities = set()
        for order in permutations(range(3)):
            if tuple(mult[i] for i in order) == mult:
                for x,y in pairs:
                    possibilities.add(tuple(sorted(tuple(v[i] for i in order) for v in (x,y))))
        assert possibilities == unshifted_pairs(mult)


def odd_part(value):
    return value // (value & -value)


def small_degrees():
    checked = removed = 0
    for alpha,beta in ((2,1),(3,1),(3,2)):
        for mult,_,_,points in CASES:
            for order in permutations(range(3)):
                numerator, zero_power = 1, 0
                for point, power in zip(points,(3,2,2)):
                    if point is None:
                        continue
                    c0,c1,c2 = (point[i] for i in order)
                    difference = c2**beta*c0**(alpha-beta)-c1**alpha
                    if difference:
                        numerator *= abs(difference)**power
                    elif point == G4:
                        removed += 1
                        break
                    else:
                        assert all(abs(x) == 1 for x in point)
                        zero_power += power
                else:
                    assert zero_power <= 3
                    assert numerator*257**zero_power < 27*256**(2*alpha)
                    checked += 1
    print('Small degree bounds:', checked, 'checked;', removed, 'discarded by the square relation')


def finite_scan():
    counts = Counter()
    for base in range(3, 256, 2):
        alpha = 2
        high = base**alpha
        while high < 2**128:
            for beta in range(1,alpha):
                if gcd(alpha,beta) != 1:
                    continue
                for values in permutations((1,base**beta,high)):
                    for mult,xf,yf,_ in CASES:
                        if base%3 == 0 and mult != (2,2,2):
                            continue
                        m,x,y = (dot(f,values) for f in (mult,xf,yf))
                        u = gcd(x,y)
                        while gcd(u,base*m) > 1:
                            u //= gcd(u,base*m)
                        s = odd_part(lcm(gcd(m,x),gcd(m,y)))
                        counts['assignments'] += 1
                        if u >= 7:
                            counts['scale at least seven'] += 1
                            assert u**3*s**2 < 27*m**2, (base,alpha,beta,values,xf,yf)
            alpha += 1
            high *= base
    print('Finite allocation search:', dict(counts))


def vector_gcd(values, vector):
    if vector is None:
        return 1
    # Actual units are in these projective ratios at the assigned odd power.
    minors = [values[i]*vector[j] - values[j]*vector[i] for i,j in ((0,1),(1,2),(2,0))]
    return odd_part(gcd(*minors))


def cross_scan():
    counts = Counter()
    for mult,_,_,points in CASES:
        for alpha in range(2, 128):
            for base in range(3,256,2):
                high = base**alpha
                if high >= 2**128:
                    break
                if base%3 == 0 and mult != (2,2,2):
                    continue
                for beta in range(1,alpha):
                    if gcd(alpha,beta) != 1:
                        continue
                    for values in permutations((1,base**beta,high)):
                        u,g1,g2 = (vector_gcd(values,v) for v in points)
                        m = dot(mult,values)
                        counts['assignments'] += 1
                        if u >= 7:
                            counts['odd bound at least seven'] += 1
                            assert u**3*(g1*g2)**2 < 27*m**2, (values,points)
    print('Independent signed-vector search:', dict(counts))


def check_repeated():
    from flint import fmpz_poly

    x = fmpz_poly([0, 1])
    counts = Counter()
    for n in range(1,29):
        modulus = x**n+1
        for a in range(2*n):
            for b in range(2*n):
                common = modulus.gcd((x**b+x**a+1) % modulus)
                t = gcd(n,a,b)
                sixth = n//t%6 == 3 and {a//t%6,b//t%6} == {2,4}
                predicted = x**(2*t)-x**t+1 if sixth else fmpz_poly([1])
                assert common == predicted
                counts['polynomial pairs'] += 1
                if common.degree():
                    counts['nonconstant gcds'] += 1
                for p in (3,5,7,11):
                    u = gcd(p**n+1,p**b+p**a+1)
                    s = odd_part(gcd(p**n+1,p**abs(b-a)+1))
                    assert u**3*s**2 < 27*p**(2*n)
                    if len({0,a%n,b%n}) == 3:
                        assert u**2 <= 3**(n-common.degree())*int(common(p))**2
                        gap = gcd(n,abs(b-a))
                        if not sixth and n == 3*gap and abs(b-a)//gap%2:
                            assert u <= 4**gap
                            counts['refined resultant checks'] += 1
    print('Repeated-part checks:', dict(counts))


def row_control():
    count = 0
    weights = set()
    for p in (3,5,7):
        for exponents in product(range(1,7),repeat=3):
            powers = tuple(p**i for i in exponents)
            for mult in ((3,1,2),(4,1,1),(2,2,2)):
                total = 1+dot(mult,powers)
                if (p,total) in weights:
                    continue
                weights.add((p,total))
                for row in free_rows(total,p):
                    count += 1
                    # For seven formal units, a scale above two requires the
                    # binary-free companion, whose normalized sum is at least seven.
                    assert gcd(*row) < 7 or has_size_obstruction(total,row)
    print('Direct allocation controls:', len(weights), 'weights,', count, 'positive rows')


if __name__ == '__main__':
    template_controls()
    small_degrees()
    finite_scan()
    cross_scan()
    check_repeated()
    row_control()
