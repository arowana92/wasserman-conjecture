"""The simultaneous local prime classes at cofactor 986, at all digit depths."""
from fractions import Fraction
from math import gcd
from check_cofactor_986_powers import ROWS, productive_edges


def spaced(v, bits):
    modulus = 1 << bits
    return ((v & (v<<1) | v & (v<<2)) % modulus) == 0


def possible_binary(w, bits):
    modulus = 1 << (bits+1)
    v = (986*w+1)*pow(7,-1,modulus) % modulus
    return spaced(v,bits+1)


def binary_classes(n):
    assert not possible_binary(1,8) and not possible_binary(-1,8)
    allowed = []
    for q in range(1,1024,2):
        if q%256 in (1,255):
            continue
        x = q**4-1
        h = (x & -x).bit_length()-1
        assert h <= 9
        if possible_binary(pow(q,n,1<<h),h):
            allowed.append(q)
    minimal = {x%128 for x in allowed}
    assert len(allowed) == 8*len(minimal) and len(minimal) == 17
    return minimal


def live_prefix(es,x,p,k):
    state = 0
    for _ in range(k):
        x,digit = divmod(x,p)
        if digit not in es[state]:
            return False
        state = es[state][digit]
    return True


def verify():
    binary = {n:binary_classes(n) for n in (1,3)}
    for (slopes, _), w, max_h in zip(ROWS,(10,21),(1,2)):
        edges = productive_edges(2,slopes)
        assert all(set(row)=={0} for row in edges.values())
        edges = productive_edges(29,slopes)
        for h in range(1,max_h+2):
            root = pow(w,29**(h-1),29**h)
            assert pow(root,28,29**h) == 1
            assert live_prefix(edges,root,29,h) == (h <= max_h)
    masses = []
    # Each primitive residue at 29 is 10**g. The congruence of the full
    # exponent at 28 then selects its binary test and its two 7-adic digits.
    for row_index,(t,max_h,targets) in enumerate(((1,1,{22,29}),(19,2,{22}))):
        mass = Fraction()
        for g in range(1,28):
            if gcd(g,28) != 1:
                continue
            n = t*pow(g,-1,28)%28
            exponent = next(k for k in range(n,n+3*28,28) if k%3 == 0)
            classes7 = {q for q in range(1,49) if q%7 and pow(q,exponent,49) in targets}
            assert len(classes7) == 3*len(targets)
            mass += (Fraction(1,28)*Fraction(29**max_h-1,29**max_h)
                     * Fraction(len(classes7),42)*Fraction(1,17)*Fraction(3,8)
                     * Fraction(len(binary[n%4]),64))
        masses.append(mass)
        print("Row",row_index+1,"relative mass",mass)
    assert masses == [Fraction(9,25984),Fraction(135,753536)]
    assert sum(masses) == Fraction(99,188384)
    for n,classes in binary.items():
        print("Exponent modulo 4:",n,"binary prime residues modulo 128:",sorted(classes))
    print("All-exponent certified relative density at least",1-sum(masses))


if __name__ == '__main__':
    verify()
