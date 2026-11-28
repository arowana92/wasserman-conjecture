"""Direct factorial oracle and an independent congruence test for affine endpoint rows."""

from itertools import combinations_with_replacement
from math import comb, gcd, prod

from check_affine_endpoint import divisors, factors, scan
from check_line_divisors import valuation


def linear_congruence(a,b,m):
    g=gcd(a,m)
    if b%g:
        return None
    a,b,m = a//g,b//g,m//g
    return b*pow(a,-1,m)%m,m


def merge(u,m,v,n):
    sol = linear_congruence(m,v-u,n)
    if sol is None:
        return None
    j,step = sol
    return (u+m*j)%(m*step),m*step


def first_by_progression(normal,e,u,s_powers,s,v):
    a,b,c = normal
    delta=b-c
    target = e//u
    if not c*(v-3)+sum(normal)<=target<=b*(v-3)+sum(normal):
        return None
    prime_bases = [p for p,q in factors(u)+s_powers]
    if a==b or a==c:
        unique = c if a==b else b
        index,rem = divmod(target-a*v,unique-a)
        if rem or not 0<index<v:
            return None
        cut=(u*index,u*(v-index))
        if any(valuation(cut,p) for p in prime_bases):
            return None
    # A power at which i cannot be the exceptional index divides part i.
    steps = [prod(q for p,q in s_powers if (e-i)%q) for i in normal]
    value,period = 0,steps[0]
    for coeff,target,modulus in ((a-c,e//u-c*v,delta*steps[1]),
                                 (a-b,e//u-b*v,delta*steps[2])):
        sol = linear_congruence(coeff,target,modulus)
        if sol is None:
            return None
        combined=merge(value,period,*sol)
        if combined is None:
            return None
        value,period=combined
    for x in range(value or period,v-1,period):
        y=(e//u-c*v-(a-c)*x)//delta
        z=v-x-y
        if min(y,z)<=0:
            continue
        row = (u*x,u*y,u*z)
        if all(valuation(row,p)==0 for p in prime_bases):
            return row
    return None


def normals():
    for a,b,c in combinations_with_replacement(range(-2,3),3):
        if a!=c and gcd(a,b,c)==1 and (a,b,c)<=(-c,-b,-a):
            yield b,c,a


def exact_oracle():
    from check_affine_endpoint import first_row
    positive = negative = 0
    for normal in normals():
        for n in range(3,27):
            for u in divisors(factors(n)):
                v=n//u
                if v<3:
                    continue
                for s in divisors(factors(n-1)):
                    possible = set()
                    for x in range(1,v-1):
                        for y in range(1,v-x):
                            row = (u*x,u*y,u*(v-x-y))
                            k = comb(n,row[0])*comb(n-row[0],row[1])
                            if gcd(k,u*s)==1:
                                possible.add(sum(a*i for a,i in zip(normal,row)))
                    for e in range(min(normal)*n,max(normal)*n+1,u):
                        expected=e in possible
                        args=(normal,e,u,factors(s),s,v)
                        assert (first_row(*args) is not None)==expected
                        assert (first_by_progression(*args) is not None)==expected
                        positive += expected
                        negative += not expected
    print(f"direct factorial decisions: positive={positive}, negative={negative}",flush=True)


if __name__=="__main__":
    exact_oracle()
    counts = 0
    for normal in normals():
        got = scan(normal,2048,10**9,row_finder=first_by_progression)
        print(normal,got['predecessor_power_test'],flush=True)
        counts += got['predecessor_power_test']
        assert not got['first_row_pass']
    print(f"progressions checked above one billion={counts}")
