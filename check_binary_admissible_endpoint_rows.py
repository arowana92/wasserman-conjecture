"""Verify an infinite two-row family passing all endpoint layers and full binary addition."""
from math import gcd, prod

from check_central_cuts import legendre
from check_scaled_row_blocks import certify_prime


C = 187
D = sum(1 << i for i in (9,47,58,89))
B = 3310000106583215246468608
FACTORS = {3:3,67:1,389:1,11071:1,424866988899517:1}


def progressions():
    assert C*B == D and gcd(C,B+1) == 1 and B%2 == 0
    assert prod(p**e for p,e in FACTORS.items()) == B-1
    for p in FACTORS:
        certify_prime(p)
        exp = (p-1)//gcd(840,p-1)
        target = -pow(C*(1 << 28),-1,p)%p
        assert pow(pow(2,840,p),exp,p) == 1
        if p == 11071:
            assert exp == 369
            assert [h for h in range(exp) if pow(2,840*h,p) == target] == [136]
        else:
            assert pow(target,exp,p) != 1
    assert pow(2,840,C*1225) == pow(2,840,9) == 1
    assert (1+C*(1 << 28))%1225 == 123 and B%1225 == 908


def concrete(h):
    assert h >= 1 and h%369 != 136
    k = 28+840*h
    z = 1 << k
    v,u = 1+C*z,1+D*z
    n,S,T = u*v,C*z,1+B*v
    l,rem = divmod(T,35)
    assert not rem and l%35 == 6 and gcd(l,35) == 1
    A = tuple(u*w for w in (18*z+1,33*z,136*z))
    row_b = tuple(v*w for w in (23*l-B,5*l,7*l))
    assert sum(A) == sum(row_b) == n and min(*A,*row_b) > 0
    assert gcd(u,v) == gcd(S,T) == 1 and S*T == n-1
    for row,m,expected in ((A,u,(0,0,0)),(row_b,v,(0,0,0)),(A,z,(1,0,0)),
                           (A,11,(0,0,1)),(A,17,(0,1,0)),(row_b,l,(1,0,0)),
                           (row_b,5,(0,0,1)),(row_b,7,(0,1,0))):
        assert tuple(x%m for x in row) == expected
    assert all((A[i] & A[j]) == 0 for i,j in ((0,1),(1,2),(0,2)))
    assert sum(bin(x).count("1") for x in A) == bin(n).count("1")
    assert legendre(n,2) == sum(legendre(x,2) for x in A)
    assert n%9 == 7 and sum(x%3 for x in A) == 4 and sum(x%9 for x in row_b) == 16
    assert all(legendre(n,3) > sum(legendre(x,3) for x in row) for row in (A,row_b))
    return n.bit_length()


if __name__ == "__main__":
    progressions()
    sizes = [concrete(h) for h in (1,2,3,4)]
    print(f"Verified the complete residue progression and weights of bit lengths {sizes}.")
