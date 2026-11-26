"""Independent finite-field controls for the two degree enumerations (standard library)."""
from collections import Counter, defaultdict
from itertools import combinations_with_replacement, permutations, product
from math import gcd, isqrt, lcm


def root_table(order):
    """A finite field containing a primitive root of this exact order."""
    p = 100000000 // order * order + 1
    while any(p % t == 0 for t in range(2, isqrt(p)+1)):
        p += order
    divisors = [t for t in range(2, order+1) if order % t == 0]
    for base in range(2, p):
        z = pow(base, (p-1)//order, p)
        if all(pow(z, order//t, p) != 1 for t in divisors):
            return p, [pow(z, k, p) for k in range(order)]
    raise AssertionError("primitive root missing")


def mixed_offsets():
    # Test oriented relative heights in a finite field by direct addition.
    options = set()
    for row in ((0,0,0,1), (0,0,0,2), (0,0,1,1), (0,0,1,2)):
        options.update(tuple(t-u[0] for t in u) for u in permutations(row))
    p, roots = root_table(12)
    radical = (roots[1]+roots[-1]) % p
    assert radical**2 % p == 3
    lut = defaultdict(set)
    for delta in options:
        weights = [pow(radical, t, p) for t in delta]
        for tail in product(range(12), repeat=3):
            phase = (0, *tail)
            values = [a*roots[k] % p for a,k in zip(weights, phase)]
            if sum(values) % p or any((a+b) % p == 0 for i,a in enumerate(values) for b in values[:i]):
                continue
            lut[phase].add(delta)
    return lut


def mixed_control():
    lut = mixed_offsets()
    counts = Counter()
    exceptions = set()
    for n in range(1, 25):
        m = lcm(2*n, 12)
        p, roots = root_table(m)
        radical = (roots[m//12]+roots[-m//12]) % p
        assert radical**2 % p == 3
        root_angles = list(range(m//(2*n), m, m//n))
        phase_columns = [[roots[k*j % m] for k in root_angles] for j in range(25)]
        for tail in combinations_with_replacement(range(25), 3):
            if gcd(n, *tail) != 1:
                continue
            slopes = (0, *tail)
            possible = set()
            for k in range(12):
                if k*n % 12 == 6:
                    possible.update(lut.get(tuple(k*j % 12 for j in slopes), ()))
            for delta in possible:
                rows = [[v*pow(radical, t, p) % p for v in phase_columns[j]]
                        for j,t in zip(slopes, delta)]
                columns = list(zip(*rows))
                sums = [sum(col) % p for col in columns]
                da = sums.count(0)
                du = [sum((s-col[i]) % p == 0 for col,s in zip(columns, sums)) for i in range(4)]
                if max(Counter(zip(slopes, delta)).values()) >= 3:
                    continue
                for b in (0, 1):
                    if b*n % 2 or any((b*j+t) % 2 for j,t in zip(slopes, delta)):
                        continue
                    if da == n:
                        counts["full"] += 1
                        continue
                    for i, degree in enumerate(du):
                        counts["allocations"] += 1
                        if 3*degree+2*da >= 2*max(n, *tail):
                            exceptions.add((n, slopes, delta, i, b, degree, da))
    assert counts == {"full": 1432, "allocations": 64008}
    assert exceptions == {
        (3, (0,1,2,4), (0,2,0,0), 1, 0, 2, 1),
        (3, (0,2,3,4), (0,0,2,0), 2, 0, 2, 1),
    }
    print("Root counts in split finite fields recover both mixed equality cases")


def polynomial(terms, p):
    powers = [j for j, _ in terms]
    a = [0]*(max(powers)-min(powers)+1)
    for j, t in terms:
        a[j-min(powers)] += pow(3, t, p)
    a = [v % p for v in a]
    assert a[0] and a[-1]
    return a


def common_degree(a, b, p):
    while b:
        a = a.copy()
        inv = pow(b[-1], -1, p)
        while len(a) >= len(b):
            c = a[-1]*inv % p
            start = len(a)-len(b)
            for k, v in enumerate(b):
                a[start+k] = (a[start+k]-c*v) % p
            while a and a[-1] == 0:
                a.pop()
        a, b = b, a
    return len(a)-1


def translation_control():
    p = 1000000007
    count = profiles = attained = 0
    for n in range(1, 10):
        for u, v in combinations_with_replacement(range(10), 2):
            if gcd(n,u,v) != 1:
                continue
            for (s,t), b in product(((0,0), (1,0), (0,1), (-1,-1)), (0,1)):
                if n*b % 2 or (u*b+s) % 2 or (v*b+t) % 2:
                    continue
                y = polynomial(((0,0), (u,(u*b+s)//2), (v,(v*b+t)//2)), p)
                f = polynomial(((0,0), (n,n*b//2)), p)
                mu = common_degree(y, f, p)
                if 100*mu < 21*max(n,v):
                    continue
                profiles += 1
                for r in range(-11, 12):
                    h = max(n, v, max(0,r+n)-min(0,r))
                    for j in range(-160, 161):
                        other = polynomial(((0,0), (r,j), (r+n,j+n*b//2)), p)
                        degree = common_degree(y, other, p)
                        assert 2*(3*mu+2*degree) <= 3*h
                        attained += 2*(3*mu+2*degree) == 3*h
                        count += 1
    assert profiles == 7 and count == 51681 and attained
    print(f"Euclidean algorithm modulo {p} confirms the bound for {count} translations")


if __name__ == "__main__":
    mixed_control()
    translation_control()
