"""Check higher-layer exclusions and the surviving perfect-power exponent degrees."""
from math import gcd, lcm

from check_binary_cofactor_rays import Shape
from check_binary_exponent_progressions import published_rays


def cost(slopes, scalar, modulus):
    x,y,z = slopes
    return ((x+y+z)*scalar+1)//modulus - (
        x*scalar//modulus + y*scalar//modulus + (z*scalar+1)//modulus)


def excluded_cofactors():
    tests = {15871:(59,3), 33661:(821,2), 62463:(443,2), 62809:(587,2)}
    for ray in published_rays():
        row = ray.row
        if row.d not in tests:
            continue
        p,depth = tests[row.d]
        first = pow(2,row.a,p)
        assert pow(2,row.h,p) == 1 and (p,p) in row.forced
        residues = [first]
        width = p
        counts = []
        for _ in range(1,depth):
            after = width*p
            residues = [t+digit*width for t in residues for digit in range(p)
                        if cost(row.slopes,t+digit*width,after) == 0]
            counts.append(len(residues))
            width = after
        assert not residues
        print(row.d, "prime",p,"scalar residue",first,"surviving lifts",counts)


def incompatible(first, period, degree, modulus):
    common = gcd(period,modulus)
    return all((pow(i,degree,modulus)-first) % common for i in range(modulus))


def local_powers(row, p, power, possible):
    """Test all compatible exponent residues at modulus p^power."""
    width = p**power
    full = lcm(row.h,width)
    assert pow(2,full,p*width) == 1
    shared = gcd(row.h,width)
    for a in possible:
        if (a-row.a) % shared:
            continue
        jump = ((a-row.a)//shared * pow(row.h//shared,-1,width//shared)) % (width//shared)
        e = row.a+row.h*jump
        assert e % width == a and e % row.h == row.a % row.h
        q = pow(2,e,p*width)
        assert any(cost(row.slopes,q,p**k) for k in range(1,power+2))


def degrees():
    for ray in published_rays():
        r = ray.row
        if r.d == 12031:
            # The refined class is derived in check_binary_exponent_progressions.py.
            a,h = 769251571,2043279108
            for degree,modulus in ((2,4),(3,13),(7,113),(13,53)):
                assert incompatible(a,h,degree,modulus)
        if r.d == 56279:
            for degree,modulus in ((2,7),(3,7),(41,83)):
                assert incompatible(r.a,r.h,degree,modulus)
            local_powers(r,7,4,{pow(i,49,7**4) for i in range(7**4)})
        if r.d == 60829:
            r = r._replace(h=627270)
            assert incompatible(r.a,r.h,2,29)
            assert r.a % 7 == r.h % 7 == 0
            local_powers(r,7,3,{0})
        if r.d == 61423:
            assert r.h % 4 == 0 and r.a % 4 == 2
        if r.d == 65095:
            assert 80500 % 25 == 0 and 40240 % 25 == 15
        if r.d == 65243:
            for degree,modulus in ((2,5),(3,13),(5,41)):
                assert incompatible(r.a,r.h,degree,modulus)
            local_powers(r,53,2,{pow(i,13,53**2) for i in range(53**2)})

    missed = ((2,3,7,13),(2,3,41,49),(2,3,5,13))
    combined = {lcm(a,b,c) for a in missed[0] for b in missed[1] for c in missed[2]}
    minimal = sorted(x for x in combined if not any(x != y and x%y == 0 for y in combined))
    assert minimal == [2,3,245,533,637,1435]
    print("Joint sufficient degrees",minimal)

    for m,q in ((2,4),(3,13),(11,56453),(1283,56453)):
        assert incompatible(2337019,2935556,m,q)
    large = Shape(7338895,3*1467779,2918400,2*1467779-2918400,3,1,2337019,2935556)
    local_powers(large,5,3,{pow(i,5,125) for i in range(125)})
    assert 3948614%4 == 2 and 5647044%4 == 0


if __name__ == "__main__":
    excluded_cofactors()
    degrees()
