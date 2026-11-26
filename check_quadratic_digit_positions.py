"""Verify six quadratic digit positions at every prime and every outer shift."""
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations, product
from math import gcd
from pathlib import Path
from random import Random
from subprocess import PIPE, Popen, run
from tempfile import TemporaryDirectory

from flint import fmpz, fmpz_poly
from check_quadrinomial_small_bases import CASES, X

INDICES = tuple(range(6))
HEIGHT_CAP = 68333
ORDER_BOUND = 100


@lru_cache(None)
def phase_tuples(p, heights):
    order = CASES[p][0][-1]
    phi = fmpz_poly.cyclotomic(order)
    phases = [(X**j) % phi for j in range(order)]
    shift = order//(8 if p == 2 else 12)
    radical = phases[shift]+phases[-shift]
    weighted = [[(z*radical**h) % phi for z in phases] for h in heights]
    answers = []
    for other in product(range(order), repeat=len(heights)-1):
        indices = (0, *other)
        if sum(block[j] for block, j in zip(weighted, indices)):
            continue
        if any(heights[i] == heights[j] and (indices[i]-indices[j]) % order == order//2
               for i, j in combinations(range(len(indices)), 2)):
            continue
        answers.append(indices)
    return answers


def phase_options(p, center, ids):
    layers = [k*(k-center) for k in ids]
    if len(ids) == 3:
        options = [(0, 0, 1), (0, 0, 2)] if p == 2 else [(0, 0, 1)]
    else:
        options = [row for row, count in CASES[p][1].items() if count]
    answers = {}
    for sign in (-1, 1):
        values = [sign*k for k in layers]
        for row in options:
            magnitude = Fraction(row[-1], 2*(max(values)-min(values)))
            heights = [2*magnitude*(k-min(values)) for k in values]
            if sorted(heights) == list(row):
                phases = phase_tuples(p, tuple(map(int, heights)))
                if phases:
                    answers[sign*magnitude] = phases
    return answers


def mixed_degree(gap, center, q, r, order, ids, phases):
    if not phases:
        return 0
    v = q+center*r
    slopes = [gap*v]+[(k-ids[0])*(q+(k+ids[0])*r) for k in ids[1:]]
    common = gcd(*slopes)
    slopes = [a//common for a in slopes]
    result = 0
    for phase in phases:
        angles = (order//2, *phase[1:])
        if all((slopes[i]*angles[j]-slopes[j]*angles[i]) % order == 0
               for i, j in combinations(range(len(angles)), 2)):
            result += common
    return result


def polynomial(p, ids, t, c, q, r):
    powers = [k*t+k*k*c for k in ids]
    slopes = [k*q+k*k*r for k in ids]
    return sum(p**(a-min(powers))*X**(b-min(slopes)) for a, b in zip(powers, slopes))


def degree_control(p, i, j, trio, q, r, determinant, du, da):
    v = q+(i+j)*r
    low = (-determinant*pow(r, -1, v)) % v if v > 1 else 0
    c0 = (determinant+r*low)//v
    t0 = low-(i+j)*c0
    assert q*c0-r*t0 == determinant
    rest = tuple(k for k in INDICES if k not in (i, j))
    binomial = polynomial(p, (i,j), t0, c0, q, r)
    assert binomial.gcd(polynomial(p, trio, t0, c0, q, r)).degree() == du
    assert binomial.gcd(polynomial(p, rest, t0, c0, q, r)).degree() == da


def degree_table():
    table = defaultdict(int)
    exceptional = set()
    tests = controls = 0
    largest_cutoff = Fraction(0)
    for i, j in combinations(INDICES, 2):
        center, gap = i+j, j-i
        rest = tuple(k for k in INDICES if k not in (i,j))
        spread = lambda ids: max(k*(k-center) for k in ids)-min(k*(k-center) for k in ids)
        error = 5*gap*spread(rest)
        margin = Fraction(2, 3)-Fraction(2*(center == 5), 5)
        ratio = max(Fraction(1, 5), Fraction(center, 25))
        assert 9*gap*ORDER_BOUND < HEIGHT_CAP*(margin-error*ratio/ORDER_BOUND)
        for p in (2, 3):
            order = CASES[p][0][-1]
            quad_options = phase_options(p, center, rest)
            for trio in combinations(rest, 3):
                tri_options = phase_options(p, center, trio)
                error = gap*(3*spread(trio)+2*spread(rest))
                for delta in sorted(set(tri_options) | set(quad_options)):
                    for v in range(1, ORDER_BOUND+1):
                        if (v*delta).denominator != 1:
                            continue
                        determinant = int(v*delta)
                        for r in range(v//center+1):
                            q = v-center*r
                            if gcd(q, r) != 1:
                                continue
                            du = mixed_degree(gap, center, q, r, order, trio, tri_options.get(delta))
                            da = v*(center == 5) + mixed_degree(
                                gap, center, q, r, order, rest, quad_options.get(delta))
                            combined = 3*du+2*da
                            key = i,j,v,r,determinant
                            if combined > 2*v*(center == 5):
                                table[key] = max(table[key], combined)
                            tests += 1
                            fraction = Fraction(combined, 5*q+25*r)
                            degree_control(p, i, j, trio, q, r, determinant, du, da)
                            controls += 1
                            if fraction >= Fraction(2, 3):
                                exceptional.add((p,i,j,trio,delta,q,r,du,da))
                                continue
                            # Multiply the upper bound for log(U^3 A^2) out
                            # using L/v = (H-(25-5*center)*delta)/(5*q+25*r).
                            constant = 9*gap*v+error*abs(determinant)-fraction*(25-5*center)*delta
                            cutoff = constant/(Fraction(2, 3)-fraction)
                            largest_cutoff = max(largest_cutoff, cutoff)
    assert largest_cutoff == Fraction(136665, 2) < HEIGHT_CAP
    bad = {(2,1,3,(0,4,5),Fraction(1,5),1,1,10,0),
           (3,1,4,(0,3,5),Fraction(-1,12),7,1,12,12),
           (3,1,4,(2,3,5),Fraction(1,12),7,1,12,12)}
    assert exceptional == bad
    for p,i,j,trio,delta,q,r,du,da in bad:
        det = int(delta*(q+(i+j)*r))
        single, = set(INDICES)-{i,j}-set(trio)
        first, second = (polynomial(p, ids, -det, 0, q, 1)
                         for ids in (trio,(i,j,single)))
        assert first.gcd(second) == 1
        assert max(first.degree(), second.degree()) <= 60
        assert max(sum(int(abs(a)) for a in f.coeffs()) for f in (first,second)) <= 3*p**5
        assert 3*du+2*da == 5*q+25*r
    assert len(table) == 72
    print(f"{tests} exceptional direction tests, {controls} independent polynomial gcd pairs, "
          f"{len(table)} finite-cover degree entries", flush=True)
    return table


def aligned_binary():
    for t, c in product(range(1, 26), range(1, 6)):
        h = [k*t+k*k*c for k in INDICES]
        for i,j,k,l,a,b in permutations(INDICES):
            if i < j and k < l and a < b and h[j]-h[i] == h[l]-h[k]:
                assert abs((h[j]-h[i])-(h[b]-h[a])) != 1


def signed_remainder(exponents, gap):
    positions = sorted((k % gap, (-1)**(k//gap)) for k in exponents)
    start = max(range(len(positions)), key=lambda j:
                (positions[j][0]-positions[j-1][0]) % gap)
    origin = positions[start][0]
    coefficients = Counter()
    for k, sign in positions:
        coefficients[(k-origin) % gap] += sign if k >= origin else -sign
    return [(e,b) for e,b in coefficients.items() if b]


def verify_integers(i, j, t, c, counts):
    h = [k*t+k*k*c for k in INDICES]
    gap, H = h[j]-h[i], h[-1]
    step = gcd(t+c, 2*c)
    if H//step <= 100:
        return
    other = [h[k] for k in INDICES if k not in (i,j)]
    quad = signed_remainder(other, gap)
    if not quad:
        return
    trios = [signed_remainder([a for a in other if a != z], gap) for z in other]
    for p in (2,3,5,7,11,13,17,19,23):
        if p**step >= 27:
            break
        counts['base comparisons'] += 1

        powers = {}

        def power(e):
            if e not in powers:
                powers[e] = fmpz(p)**e
            return powers[e]

        def divisor(terms):
            if len(terms) == 1:
                return fmpz(abs(terms[0][1]))
            return (power(gap)+1).gcd(sum(b*power(e) for e,b in terms))

        A = divisor(quad)
        if len(quad) > 1 and A == power(gap)+1:
            continue
        for z, terms in zip(other, trios):
            if terms and 9*(max(e for e,b in terms)+2)+6*A.bit_length() <= 2*H:
                continue
            U = divisor(terms) if terms else power(gap)+1
            if 9*U.bit_length()+6*A.bit_length() <= 2*H:
                continue
            if 16**3*U**9*A**6 < 27**3*power(2*H):
                continue
            y = sum(power(e) for e in other if e != z)
            x = power(h[i])+power(h[j])
            C = y.gcd(x+power(z))
            assert U**3*A**2*C**2 < 27*(x+y+power(z))**2, (p,i,j,t,c,z)
            counts['full product comparisons'] += 1



def direct_product_controls(i, j, t, c):
    # Independently evaluate the original sums on a deterministic sample.
    # Every C divisor is coprime to p, so shifted monomials can be removed.
    h = [k*t+k*k*c for k in INDICES]
    step = gcd(t+c, 2*c)
    rest = tuple(k for k in INDICES if k not in (i,j))
    count = 0
    for p in (2,3,5,7,11,13,17,19,23):
        if p**step >= 27 or h[-1]//step <= 100:
            break
        values = [fmpz(p)**a for a in h]
        total = sum(values)
        f = 1+fmpz(p)**(h[j]-h[i])
        a = f.gcd(sum(values[k] for k in rest))
        for single in rest:
            trio = [k for k in rest if k != single]
            y = sum(values[k] for k in trio)
            u = f.gcd(y)
            if u < (3 if p == 2 else 7):
                continue
            other = (i,j,single)
            c_bound = (y//min(values[k] for k in trio)).gcd(
                sum(values[k] for k in other)//min(values[k] for k in other))
            assert u**3*a**2*c_bound**2 < 27*total**2
            count += 1
    return count

def sampled_cover(table, i, j, t, c):
    if gcd(t,c) >= 5:
        return True
    B, n, H = i+j, j-i, 5*t+25*c
    heights = [k*(k-B) for k in INDICES if k not in (i,j)]
    error = 5*n*(max(heights)-min(heights))
    for v in range(1, ORDER_BOUND+1):
        for r in range(v//B+1):
            q = v-B*r
            delta = c*q-t*r
            if not delta or gcd(q,r) != 1:
                continue
            degree = table.get((i,j,v,r,delta), 2*v*(B == 5))
            if 3*(degree*(t+B*c)+(9*n*v+error*abs(delta))*v) <= 2*H*v:
                return True
    return False


def finite_check(table):
    rng = Random(785)
    sample = set()
    for i,j in combinations(INDICES, 2):
        for _ in range(80):
            c = rng.randint(1, (HEIGHT_CAP-5)//25)
            t = rng.randint(1, (HEIGHT_CAP-25*c)//5)
            sample.add((i,j,t,c))
    root = Path(__file__).resolve().parent
    counts = Counter()
    with TemporaryDirectory(dir=root) as tmp:
        exe = Path(tmp)/"cover"
        run(["c++","-O3","-std=c++17",str(root/"check_quadratic_digit_cover.cpp"),
             "-o",str(exe)], check=True)
        with Popen([str(exe)], stdin=PIPE, stdout=PIPE, text=True) as process:
            for key, value in sorted(table.items()):
                print(*key, value, file=process.stdin)
            process.stdin.close()
            for line in process.stdout:
                case = tuple(map(int, line.split()))
                sample.discard(case)
                verify_integers(*case, counts)
                counts['residual allocations'] += 1
                if counts['residual allocations'] % 199 == 0:
                    counts['independent full products'] += direct_product_controls(*case)
            assert process.wait() == 0
    for case in sample:
        assert sampled_cover(table, *case), case
    assert counts['residual allocations'] == 411088
    print(dict(counts), f"{len(sample)} independently recertified interval exclusions")


if __name__ == '__main__':
    aligned_binary()
    finite_check(degree_table())
