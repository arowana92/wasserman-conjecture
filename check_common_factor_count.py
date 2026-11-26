"""Verify the parameter reconstruction, entropy margins, and local CRT counts."""
from fractions import Fraction as F
from itertools import product
from math import comb, gcd, prod

from check_binary_windows import LN2_HI, log2_upper, tangent
from search_endpoint import atoms, smallest_factors


def exponent_checks():
    rows = (
        (3, 3, 80910, 70108, 63946, 71319, 89460),
        (2, 3, 85117, 55158, 59869, 62045, 96200),
        (3, 2, 75132, 67204, 63107, 69422, 91900),
        (2, 2, 78673, 46489, 57935, 59279, 97800),
    )
    for a,b,*values in rows:
        k, lam, r0, s0, target = (F(v,100000) for v in values)
        if b == 3:
            small = (8*k-2)/5
        else:
            small = (5*k-1)/3
        assert max(small,2*(1-k)) < target
        for r,s in product((0,1), repeat=2):
            edge = (2-lam*(1+F(1,a)))*tangent(r,r0)+2*lam*tangent(F(r,a),r0/a)
            center = (1-lam)*tangent(s,s0)+lam*tangent(F(s,a),s0/a)
            assert (1-k)*edge+(2*k-1)*center < target
    theta, k, delta = F(37,50), F(7124,10000), F(1,100000)
    b = tangent(F(1,3), F(1,3))
    cost = ((1-theta)*b+theta*(tangent(F(1,6), F(1,6))+log2_upper(F(5))/1200)
            + F(2,3)*(theta-k))
    assert max(cost,(8*k-2)/5) < theta-delta
    assert 2*theta/(400**2*LN2_HI) > delta
    print("Four global entropy margins and the two-ended interval margins passed")


def check_tuple(u,v,s,d,t,z,y):
    n = u*v
    assert n == 1+s*d*t and u == z*t+d*y
    lifted = 1+s*d*pow(v*z,-1,d) if d > 1 else 1
    assert n % (s*d*d) == lifted % (s*d*d)
    modulus = v*s*d*d
    shift = (-lifted * pow(s*d*d,-1,v)) % v
    residue = lifted+shift*s*d*d
    assert n % modulus == residue % modulus
    assert s*d*d*y-z == (s*d-v*z)*u > 0
    assert (u-d*y)//z == t
    return n % modulus


def row_checks():
    # Start with actual positive rows, and infer the bins from their residues.
    limit = 55
    sieve = smallest_factors(limit**2)
    count = 0
    for u in range(3,limit+1):
        for v in range(2,limit+1):
            if gcd(u,v) > 1:
                continue
            powers = [w for _,w in atoms(u*v-1,sieve)]
            for a in range(1,u//3+1):
                for b in range(a,(u-a)//2+1):
                    row = (a,b,u-a-b)
                    assigned = [w for w in powers if sum(x%w==0 for x in row)==2]
                    total = prod(assigned)
                    d,x = min((total//gcd(x,total),x) for x in row)
                    t = total//d
                    assert x % t == 0 and (u-x) % d == 0
                    check_tuple(u,v,(u*v-1)//total,d,t,x//t,(u-x)//d)
                    count += 1
    # Generate a progression directly so the lifted classes have large periods
    # and multiple distinct weights occupy a single parameter class.
    progressions = 0
    for s,v,d,z in product(range(1,9),range(2,10),range(2,14),range(1,9)):
        if gcd(v,s*d)>1 or gcd(z,d)>1 or s*d <= v*z:
            continue
        t0 = pow(v*z,-1,d)
        t0 += d*(((-pow(s*d,-1,v)-t0)*pow(d,-1,v)) % v)
        last = None
        for step in (1,2,3):
            t = t0+step*v*d
            u = (s*d*t+1)//v
            residue = check_tuple(u,v,s,d,t,z,(u-z*t)//d)
            assert last is None or last == residue
            last = residue
        progressions += 1
    n = 5752
    left = comb(n,568)
    right = comb(n,719)*comb(n-719,719)
    assert gcd(n*(n-1),gcd(left,right)) == 1 and left%2 == 1
    print(f"{count:,} complementary-row assignments, {progressions:,} lifted progressions, and the exact mixed witness")


def cylinder_checks():
    tested = 0
    for L in range(5,11):
        for b in range(L):
            for length in range(1,L-b+1):
                for k in range(2,25):
                    groups = {}
                    e = (k & -k).bit_length()-1
                    for n in range(1<<(L-1),1<<L):
                        key = (n>>(b+length),n % (1<<b),n % k)
                        groups[key] = groups.get(key,0)+1
                    period = k//(1<<min(e,b))
                    assert max(groups.values()) <= ((1<<length)-1)//period+1
                    tested += 1
    print(f"{tested:,} finite interval/modulus tests with both ends prescribed")


if __name__ == "__main__":
    exponent_checks()
    row_checks()
    cylinder_checks()
