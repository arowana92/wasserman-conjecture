"""Independent finite checks of the rational-line theorems (standard library)."""

from collections import Counter, deque
from itertools import product
from math import comb, factorial, gcd

from check_line_divisors import line_gcd


def digit_sum(n, p):
    value = 0
    while n:
        n, digit = divmod(n, p)
        value += digit
    return value


def factorial_valuation(n, p):
    value = 0
    while n >= p:
        n //= p
        value += n
    return value


def row_val(row, p):
    return factorial_valuation(sum(row), p) - sum(factorial_valuation(x,p) for x in row)


def folded_classes():
    classes = rows = attained = 0
    for p, max_h in [(2,5), (3,3), (5,2), (7,2)]:
        for h, k in product(range(1,max_h+1), repeat=2):
            d = gcd(h,k)
            m = p**(h+k)-1
            divisor = p**d-1
            a, b = (p**h-1)//divisor, (p**k-1)//divisor
            frequency = Counter()
            for middle in range(divisor,m+1,divisor):
                u = middle//divisor
                t = (-pow(a,-1,b)*u)%b or b
                z = (u+a*t)//b
                row = (p**h*t,middle,z)
                weight = sum(row)
                assert (p**h-1)*row[0]+p**h*row[1] == p**h*(p**k-1)*row[2]
                assert (-(p**h-1)*weight)%m == middle%m
                value = row_val(row,p)
                assert value*(p-1) == digit_sum(middle,p)
                frequency[value] += 1
                attained += 1
            assert min(frequency)==d and max(frequency)==h+k and frequency[h+k]==1
            assert frequency[d] == comb((h+k)//d+p-2,p-1)**d
            classes += m//divisor
            for t, z in product(range(1,36),repeat=2):
                middle = (p**k-1)*z-(p**h-1)*t
                if middle<=0:
                    continue
                row = (p**h*t,middle,z)
                weight = sum(row)
                residue = (-(p**h-1)*weight-1)%m+1
                val = row_val(row,p)
                assert val*(p-1) >= digit_sum(residue,p)
                assert val == digit_sum(middle,p)//(p-1)+row_val((t,p**k*z),p)
                rows += 1
    print(f"sharp classes={classes}, constructed minima={attained}, sampled rows={rows}")


def parallel_offsets():
    classes = shifted = 0
    for p, bound in [(2,4),(3,3),(5,2)]:
        for h,k in product(range(1,bound+1),repeat=2):
            d = gcd(h,k)
            divisor = p**d-1
            m = p**(h+k)-1
            a, b = (p**h-1)//divisor, (p**k-1)//divisor
            for u in range(p**h):
                for n_class in range(m//divisor):
                    y = (-(p**h-1)*n_class-u-1)%m+1
                    t = (-pow(a,-1,b)*((y+u)//divisor))%b or b
                    z = ((y+u)//divisor + a*t)//b
                    row = (p**h*t+u,y,z)
                    n = sum(row)
                    bound_val, rem = divmod(digit_sum(u,p)+digit_sum(y,p),p-1)
                    assert rem==0 and bound_val>=d
                    assert n%(m//divisor)==n_class
                    assert p**h*n-row[0]-p**(h+k)*z==-u
                    assert min(row)>0 and row_val(row,p)==bound_val
                    block = 1
                    while block<=n:
                        block *= p
                    larger = (row[0]+p**h*(p**k-1)*block,y,z+(p**h-1)*block)
                    assert sum(larger)-n==m*block
                    assert p**h*sum(larger)-larger[0]-p**(h+k)*larger[2]==-u
                    assert row_val(larger,p)==bound_val
                    classes += 1
                    shifted += 1
    print(f"offset/class minima={classes}, disjoint high-block checks={shifted}")


def no_carry_row(r,s,t,p):
    # Deliberately use all three coordinates, with no edge reductions.
    transitions = [[] for _ in range(p)]
    for i in range(p):
        for j in range(p-i):
            for k in range(p-i-j):
                shift = r*i+s*j-t*k
                mask = int(i>0)+2*(j>0)+4*(k>0)
                transitions[(-shift)%p].append((shift,mask))
    todo = deque([(0,0)])
    seen = set(todo)
    while todo:
        carry, mask = todo.popleft()
        for shift, occupied in transitions[carry%p]:
            state = ((carry+shift)//p,mask|occupied)
            if state==(0,7):
                return True
            if state not in seen:
                seen.add(state)
                todo.append(state)
    return False


def power(n,p):
    if n<=1:
        return False
    while n%p==0:
        n//=p
    return n==1


def check_census():
    cases = 0
    for r, s, t in product(range(1,25),repeat=3):
        if r>s or gcd(gcd(r,s),t)>1:
            continue
        for p in [2,3,5,7]:
            forced = gcd(r+t,s+t)%p==0 or (
                s==r+1 and power(s,p) and t%s==0 and power(1+t//s,p))
            assert no_carry_row(r,s,t,p) != forced, (r,s,t,p)
            cases += 1
    print(f"full-mask carry-free decisions={cases}")


def check_gcd():
    cases = exact = 0
    f = [factorial(n) for n in range(201)]
    for r,s,t in product(range(1,8),repeat=3):
        if r>s or gcd(gcd(r,s),t)>1:
            continue
        divisor, witnesses = line_gcd(r,s,t)
        running = 0
        for x in range(1,120):
            for y in range(1,120):
                z, rem = divmod(r*x+s*y,t)
                n = x+y+z
                if rem or n>200:
                    continue
                coefficient = f[n]//(f[x]*f[y]*f[z])
                assert coefficient%divisor==0
                running = gcd(running,coefficient)
                exact += 1
        for p,(exponent,row) in witnesses.items():
            assert row_val(row,p)==exponent and min(row)>0
            assert r*row[0]+s*row[1]==t*row[2]
        # Continue until a finite list of exact coefficients attains the gcd.
        # Enumerate by weight, solving for the second coordinate directly.
        for n in range(3,512):
            if running==divisor:
                break
            for x in range(1,n-1):
                y, rem = divmod(t*n-(t+r)*x,s+t)
                if rem or not 0<y<n-x:
                    continue
                running = gcd(running,comb(n,x)*comb(n-x,y))
                exact += 1
                if running==divisor:
                    break
        assert running==divisor, (r,s,t,running,divisor)
        cases += 1
    print(f"line gcds from direct factorials={cases}, exact coefficients={exact}")


if __name__ == "__main__":
    folded_classes()
    parallel_offsets()
    check_census()
    check_gcd()
