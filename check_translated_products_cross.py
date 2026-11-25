"""Check root interpolation, carry obstructions, the adjacent-gap test, and finite row counts."""

from itertools import combinations, product
from math import comb, gcd, prod

from sympy import Poly, symbols
from search_endpoint import smallest_factors

X = symbols("X")


def shifted_polynomials(hs):
    F = Poly(prod(X+h for h in hs),X)
    basis = []
    for h in hs:
        quotient = F.exquo(Poly(X+h,X))
        basis.append(quotient.to_field().mul_ground(1/quotient.eval(-h)))
    return F,basis


def digits_carry(parts, p):
    """Use the factorial floors at each layer independently."""
    n = sum(parts)
    w = p
    while w <= n:
        if n//w > sum(x//w for x in parts):
            return True
        w *= p
    return False


def sign_checks():
    root_vectors = mixed = shifted = unshifted = 0
    for hs in ((0,1),(0,2),(0,1,2),(0,1,3),(0,2,5),
               (0,1,2,3),(0,1,3,4,6),(0,1,2,3,4,5)):
        n = len(hs)
        F,basis = shifted_polynomials(hs)
        factorizations = [Poly(d*F+eps,X).factor_list()[1]
                          for d in range(1,11) for eps in (-1,1)]
        assert all(len(factors) == 1 for factors in factorizations)
        for i, ((_,e),) in enumerate(factorizations):
            assert e in ((1,) if i % 2 == 0 else (1,2))
        for assignment in product(range(3),repeat=n):
            remainders = [sum((b for b,k in zip(basis,assignment) if k == i),Poly(0,X))
                          for i in range(3)]
            if any(c.q != 1 for f in remainders for c in f.all_coeffs()):
                continue
            root_vectors += 1
            if len(set(assignment)) > 1:
                mixed += 1
                assert hs in ((0,1),(0,1,2))
            for i,f in enumerate(remainders):
                assert all(f.eval(-h) == int(k==i) for h,k in zip(hs,assignment))
            # Give each part a positive leading term. Vary its scale so
            # equal parts or repeated low digits are not a necessary reason for rejection.
            parts = [v*F+r for v,r in zip((1,2,4),remainders)]
            shifted += check_layers(parts,hs[-1])
        for size in range(1,n+1):
            if 3*((size)//2)+1+(size-1)//2 <= n:
                continue
            for selected in combinations(hs,size):
                H = Poly(prod(X+h for h in selected),X)
                parts = [H,2*H,7*F-3*H]
                unshifted += check_layers(parts,selected[-1])
    assert mixed and not digits_carry((5,10,5),5)
    assert comb(20,5)*comb(15,10) % 5 != 0
    for n in range(2,301):
        maximum = max(b for b in range(1,n+1) if 3*(b//2)+1+(b-1)//2 <= n)
        assert maximum == n//2+int(n%4==1)
        assert 3*(2*n*maximum+n-maximum) < 4*n*n
    print(f"{root_vectors} integral root patterns ({mixed} mixed), {shifted} shifted and "
          f"{unshifted} unshifted carry checks; all degree margins through 300 agree.")


def check_layers(parts,h):
    translated = [Poly(f.as_expr().subs(X,X-h),X) for f in parts]
    n = max(f.degree() for f in parts)
    count = 0
    for t in (1,3):
        C = max(abs(int(f.nth(j)))*t**j for f in translated+[sum(translated,Poly(0,X))]
                for j in range(n))
        w = 5
        while w <= max(h+2,8*(C+1)):
            w *= 5
        entries = [int(f.eval(t*w-h)) for f in parts]
        assert min(entries) > 0 and digits_carry(entries,w) and digits_carry(entries,5)
        count += 1
    return count


def gap_checks(bound=300):
    spf = smallest_factors(bound)
    primes = [(p,p) for p in range(2,bound+1) if spf[p] == p]
    powers = []
    for p,w in primes:
        while w <= bound:
            powers.append((w,p))
            w *= p
    counts = [0,0,0]
    for ((x,p),(y,q),(z,r)) in combinations(sorted(powers),3):
        if len({p,q,r}) < 3:
            continue
        M = x*y*z
        parts = [(M//w)*pow(M//w,-1,w) for w in (x,y,z)]
        if sum(parts) != M+1 or gcd(*parts) < 3:
            continue
        counts[0] += 1
        h,k = y-x,z-y
        assert (sum(parts)-max(parts))*h*k*x >= M*(x-1)
        if x > h*k*(2*h+k):
            assert digits_carry(parts,p)
            counts[1] += 1
        elif not digits_carry(parts,p):
            counts[2] += 1
    assert all(counts)
    print(f"CRT cross-check: {counts[0]} eligible rows; {counts[1]} in the excluded gap region;"
          f" {counts[2]} no-carry controls elsewhere.")


def count_shapes():
    total = 0
    for shifts in (*combinations(range(1,6),2),(1,2,3),(1,2,4)):
        G = Poly(prod(X+h for h in shifts),X)
        digits = list(map(int,G.all_coeffs()[1:]))
        ordered = prod(comb(k+2,2) for k in digits)-2*prod(k+1 for k in digits)+1
        fixed = prod(k//2+1 for k in digits)-1
        # Burnside for swapping the two slots without the top monic unit.
        total += (3*ordered+fixed)//2
        if len(shifts) == 2:
            total += (ordered+fixed)//2
    assert total == 430530
    print("Independent coefficient-allocation count:",total)


if __name__ == "__main__":
    sign_checks()
    gap_checks()
    count_shapes()
