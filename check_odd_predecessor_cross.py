"""Independent exponent coverage using factorial floors and complete scalar prefix trees."""

from argparse import ArgumentParser
from collections import Counter, deque
from functools import lru_cache
from math import gcd, lcm
from pathlib import Path

from check_three_factors import powers


def table(path, headers):
    return {tuple(map(int, ln.split())) for ln in path.read_text().splitlines()
            if ln and (not headers or ln[0].isdigit())}


@lru_cache(maxsize=32)
def inverse_orbit(mod):
    residue = 1
    orbit = {}
    while residue not in orbit:
        orbit[residue] = len(orbit)
        residue = (residue + (residue & 1)*mod)//2
    return orbit


def intersection(i, m, j, n):
    period = lcm(m, n)
    if (i-j) % gcd(m, n):
        return None
    order = period//m
    reduced_mod = n//order
    offset = (j-i)//reduced_mod*pow(m//reduced_mod, -1, order)
    return (i+m*offset) % period, period


def lines_from_shapes(shapes):
    rows = set()
    for p, q, a, b in sorted(shapes, key=lambda v: v[1]):
        x, y = p*a, q*b
        ip, iq = inverse_orbit(p), inverse_orbit(q)
        if x%q not in iq:
            continue
        pair = intersection(ip[y%p], len(ip), iq[x%q], len(iq))
        if pair is None:
            continue
        k, step = pair
        # Compare individual actual parts to the whole weight, modulo the square.
        if any((s*pow(2, k, f*f))%(f*f) > (1+p*q*pow(2,k,f*f))%(f*f)
               for s, f in ((x,p), (y,q))):
            continue
        for r, top in powers(gcd(a,b)).items():
            if r == 2:
                continue
            w = r
            while w <= top:
                orbit = inverse_orbit(w)
                target = -p*q % w
                if target in orbit:
                    intersect = intersection(k, step, orbit[target], len(orbit))
                    if intersect:
                        u, H = intersect
                        scale = pow(2, u, w*r)
                        if ((x*scale)%(w*r) + (y*scale)%(w*r)) < w*r:
                            rows.add((p,q,a,b,r,w,u,H))
                w *= r
    return rows


def floors_ok(slopes, P, w):
    x, y, z = slopes
    return ((sum(slopes)*P+1)//w ==
            (x*P)//w + (y*P)//w + (z*P+1)//w)


def allowed(slopes, P, modulus, f, assigned):
    if assigned and modulus > assigned and (sum(slopes)*P+1) % (assigned*f) == 0:
        return False
    while modulus > 1:
        if not floors_ok(slopes, P, modulus):
            return False
        modulus //= f
    return True


def exclude_infinite(row, f, max_nodes):
    p,q,a,b,r,w,k,H = row
    slopes = (p*a,q*b,p*q-p*a-q*b)
    # Find all possible next exponents by their actual order at a prefix.
    mod = f
    while pow(2, H, mod*f) == 1 or (f == r and mod <= w):
        mod *= f
    prefixes = []
    t = pow(2, k, mod)
    step = pow(2, H, mod)
    initial = t
    while True:
        if allowed(slopes, t, mod, f, w if f == r else 0):
            prefixes.append(t)
        t = t*step % mod
        if t == initial:
            break
    # Equal triples of integer quotients have exactly the same possible suffixes.
    # Explore representatives as actual scalar prefixes. Legendre floors avoid
    # invoking the digit-by-digit transition formula of the primary verifier.
    todo = deque((v,mod) for v in prefixes)
    visited = set()
    while todo:
        v, modulus = todo.popleft()
        state = tuple((s*v+delta)//modulus for s,delta in zip(slopes,(0,0,1)))
        if state == (0,0,0):
            return False
        if state in visited:
            continue
        visited.add(state)
        if len(visited) >= max_nodes or f > max_nodes:
            return False
        upper = f*modulus
        for t in range(v, upper, modulus):
            if floors_ok(slopes, t, upper):
                todo.append((t,upper))
    return True


def bounded(row, limit, counts):
    p,q,a,b,r,w,k,H = row
    slopes = (p*a,q*b,p*q-p*a-q*b)
    # Split the exponent interval only until a modest number of individual
    # exponents remain. Each number is then checked using modular factorial floors.
    stack = [(k or H,H,r)]
    while stack:
        start, step, modulus = stack.pop()
        if start > limit:
            continue
        if not allowed(slopes, pow(2, start, modulus), modulus, r, w):
            continue
        if (limit-start)//step < 8:
            for exponent in range(start,limit+1,step):
                counts["individual exponents"] += 1
                for depth in (4,8,16,32,64,128):
                    if any(not allowed(slopes, pow(2,exponent,f**depth), f**depth,
                                       f,w if f==r else 0) for f in (r,p,q)):
                        break
                else:
                    raise ArithmeticError(("unresolved exponent", row, exponent))
        else:
            upper = modulus*r
            period = step
            if pow(2, period, upper) != 1:
                period *= r
            assert pow(2, period, upper) == 1
            for n in range(start,min(start+period,limit+1),step):
                stack.append((n,period,upper))


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("shapes",type=Path)
    parser.add_argument("lines",type=Path)
    parser.add_argument("survivors",type=Path)
    parser.add_argument("exclusions",type=Path)
    parser.add_argument("--bound",type=int,default=10**12)
    args = parser.parse_args()
    rows = table(args.lines,True)
    assert lines_from_shapes(table(args.shapes,False)) == rows
    survivors = table(args.survivors,True)
    rejected = {s[:-1]:s[-1] for s in table(args.exclusions,False)}
    assert survivors <= rows and set(rejected) == rows-survivors
    counts = Counter()
    for row in sorted(rows):
        p,q,_,_,r,_,_,_ = row
        if row not in survivors:
            f = rejected[row]
            assert f in (r,p,q) and exclude_infinite(row,f,15000), row
            counts["unbounded exclusions"] += 1
        else:
            bounded(row,args.bound,counts)
            counts["bounded progressions"] += 1
    # For a prime exponent sharing a factor with H, that prime is at most H.
    prime_lines = [row for row in survivors if gcd(*row[-2:]) == 1]
    assert all(0 < row[-2] < row[-1] <= args.bound for row in survivors)
    print("prime exponent classes", len(prime_lines), "pairs", len({r[:2] for r in prime_lines}))
    for *_, e, H in survivors:
        residue, capped, j = 1, 1, 1
        while residue:
            if capped > args.bound and residue == e:
                raise ArithmeticError(("factorial exponent", j, e, H))
            j += 1
            residue = (residue*j) % H
            capped = min(capped*j, args.bound+1)
    print(dict(counts))
