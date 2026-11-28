"""Independent checks of uniform bounds for affine-line gcds."""

from itertools import product
from math import comb, gcd

from check_line_divisors import line_gcd, minimum, prime_factors, valuation


def cut_gcd(t,e):
    if e==-1:
        return t
    if e==0:
        return 4 if t==1 else t+1
    if e<0:
        return 1
    primes = prime_factors(t+1)
    if len(primes)!=1:
        return 1
    p, = primes
    while e%p==0:
        e//=p
    return p if e==1 else 1


def cut_tests():
    count = exact = 0
    for t,e in product(range(1,17),range(-40,41)):
        g = cut_gcd(t,e)
        assert line_gcd(1,1,t,e)[0]==g,(t,e)
        sample = 0
        for z in range(1,241):
            n = (t+1)*z+e
            if n-z<2:
                continue
            # The gcd of the interior binomial coefficients in row n-z is
            # also checked here by actual coefficients, stopping at equality.
            prefix = comb(n,z)
            for x in range(1,(n-z)//2+1):
                coefficient = prefix*comb(n-z,x)
                assert coefficient%g==0
                sample = gcd(sample,coefficient)
                exact += 1
                if sample==g:
                    break
            if sample==g:
                break
        assert sample==g,(t,e,sample,g)
        count += 1
    print(f"cut-line gcds={count}, direct multinomials={exact}")


def support_tests():
    checks = direct = 0
    for r,s,t in product(range(1,6),repeat=3):
        if r>s or gcd(r,s,t)!=1:
            continue
        for p in (2,3,5,7,11,13):
            axis = min((v for v in (r,s,-t) if v%p),key=abs)
            endpoints = sorted((0,-axis))
            max_val = max(minimum(r,s,t,p,e)[0]
                          for e in range(endpoints[0],endpoints[1]+1))
            assert p<=s+t or max_val==0
            for e in (-89,-37,-12,8,34,83):
                assert minimum(r,s,t,p,e)[0]<=max_val
                checks += 1
            direct += 1
    print(f"prime/core windows={direct}, out-of-window valuations={checks}")


def templates(v,p,modulus):
    # Find positive carry-free rows covering every residue. Multipliers chosen
    # congruent to one keep residues intact while separating digit blocks.
    step = p
    while step%modulus != 1%modulus:
        step *= p
    shift = step
    while shift<=modulus:
        shift *= step
    found = {}
    for x,y,z in product(range(1,modulus+1),repeat=3):
        row = (x,shift*y,shift**2*z)
        f = sum(a*b for a,b in zip(v,row))
        found.setdefault(f%modulus,(f,row))
    assert len(found)==modulus
    return found


def digit_block_tests():
    checks = alternatives = 0
    for v in [(1,2,-6),(4,6,-5),(2,3,-4),(3,5,-6),(6,5,-3)]:
        for p,idx in product((2,3,5),(0,2)):
            a = v[idx]
            bidx = next(i for i,x in enumerate(v) if x%p)
            b = v[bidx]
            m = abs(a)
            while m%p==0:
                m//=p
            qpart = abs(a)//m
            by_class = templates(v,p,m)
            size = 1
            while size<=max(sum(row) for f,row in by_class.values()):
                size *= p
            base = size*qpart
            for low_power in [1,p,p**3,p**11]:
                for low in {0,low_power//2,low_power-1}:
                    for residue in range(m):
                        f,row = by_class[residue]
                        # Two different high continuations of the selected word.
                        for step in (size,base*p):
                            remaining = f+a*step
                            offset = b*low+low_power*remaining
                            word = (pow(b,-1,base)*f)%base
                            digit_block = (pow(b,-1,low_power*base)*offset)%(low_power*base)
                            assert digit_block==low+low_power*word
                            assert remaining%m==residue
                            whole = [x*low_power for x in row]
                            whole[idx] += low_power*step
                            whole[bidx] += low
                            assert min(whole)>0 and valuation(whole,p)==0
                            assert sum(x*y for x,y in zip(v,whole))==offset
                            checks += 1
                        alternatives += 1
    print(f"selected digit blocks={alternatives}, positive carry-free witnesses={checks}")


if __name__=="__main__":
    cut_tests()
    support_tests()
    digit_block_tests()
