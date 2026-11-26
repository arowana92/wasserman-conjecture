"""Check the general two-power test, quartic-degree cuts, and a five-factor family."""
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from itertools import combinations
from math import comb, gcd, prod

from check_cluster_certificates import is_prime, universal_prime
from check_cluster_minima import candidate_bases
from check_quartic_triangle import binomial_order, common_cuts, digits
from check_scaled_row_blocks import certify_prime
from search_endpoint import atoms, smallest_factors


def signed_digits(value, base):
    result = []
    while value:
        value, r = divmod(value, base)
        if 2*r > base:
            value += 1
            r -= base
        result.append(r)
    return result


def pair_test(n, w, z):
    assert w < z and gcd(w, z) == 1 and n % (w*z) == 0
    coef = digits(n//w, w)
    degree = len(coef)-1
    u = z-w
    bound = degree*(1+u+sum(c*(1+u)**j for j, c in enumerate(coef) if j))
    if 2*bound+2 >= z:
        return False
    target = signed_digits(n//z, z)
    target += [0]*(degree+1-len(target))
    assert sum(t*z**j for j, t in enumerate(target)) == n//z
    return (len(target) == degree+1 and all(2*abs(t)+2 < z for t in target)
            and sum(t <= 0 for t in target) <= 1)


def part_capacity(n, w, z):
    """Upper bound on the number of positive parts jointly avoiding both bases."""
    assert w < z and gcd(w, z) == 1 and n % (w*z) == 0
    d = digits(n//w, w)
    e = len(d)-1
    u = z-w
    correction_lo = -(sum(d[j]*u**j for j in range(1, e+1, 2))//z)
    correction_hi = (d[0]+sum(d[j]*u**j for j in range(2, e+1, 2)))//z
    ranges = [(correction_lo, d[1]+correction_hi)] + [(0, t) for t in d[2:]]
    target = digits(n//z, z)
    capacity = 0
    for k in range(e+1):
        lo = hi = 0
        for j in range(max(1, k), e+1):
            entry = (-u)**(j-k)*(comb(j+1, k+1) if k else j)
            x, y = ranges[j-1]
            lo += min(entry*x, entry*y)
            hi += max(entry*x, entry*y)
        if hi >= z or lo <= 1-z:
            return None
        if lo < 0 and k < len(target):
            capacity += target[k]//(z+lo-1)
    return capacity


def small_controls():
    primes = smallest_factors(200)
    powers = [(w, factors[0][0]) for w in range(2, 201)
              if len(factors := atoms(w, primes)) == 1]
    pairs = bins = base_powers = improved = 0
    for (w, p), (z, q) in combinations(powers, 2):
        if z-w > 20 or p == q:
            continue
        for k in range(2, 101):
            n = w*z*k
            simple = pair_test(n, w, z)
            bound = part_capacity(n, w, z)
            if not simple and (bound is None or bound >= 2):
                continue
            pairs += 1
            base_powers += p != w or q != z
            improved += not simple
            for cut in range(w*z, n, w*z):
                assert binomial_order(n, cut, p) or binomial_order(n, cut, q)
                bins += 1
    assert base_powers and bins
    assert improved
    print(f"{pairs} accepted tests, {base_powers} using a higher power, {bins} exact cut checks, "
          f"{improved} exclusions by the position bounds alone.")


def maximum_parts_controls():
    powers = [(w, base) for base in (2, 3, 5, 7, 11, 13, 17, 19)
              for w in (base, base**2) if w < 150]
    checked = attained = multiple = 0
    for (w, p), (z, q) in combinations(sorted(powers), 2):
        if p == q:
            continue
        block = w*z
        longest = [0]
        for k in range(1, 161):
            n = k*block
            best = 1
            # Concatenate two optimal rows whenever their binomial split
            # avoids both primes. Conversely every split of an avoiding row
            # avoids both, and each part is a multiple of w*z.
            for i in range(1, k//2+1):
                cut = i*block
                if binomial_order(n, cut, p) == binomial_order(n, cut, q) == 0:
                    best = max(best, longest[i]+longest[k-i])
            longest.append(best)
            cap = part_capacity(n, w, z)
            if cap is not None:
                assert best <= cap
                checked += 1
                attained += best == cap
                multiple += best > 1
    assert attained and multiple
    print(f"{checked} exact maximum-length comparisons, {attained} attained bounds, "
          f"{multiple} maxima above one.")


def quartic_degree():
    w, z = 65536, 65537
    n = w*z*(z+4)**3
    certify_prime(z)
    assert pair_test(n, w, z)
    tested, avoiding = common_cuts(n, w, z, 2, z)
    assert tested and not avoiding
    print(f"Degree four at 65536, 65537: all {tested} surviving proper cuts excluded.")


def five_factors():
    p = 350708197319
    gaps = (0, 2, 12, 72, 102)
    a, b = gaps[1:3]
    s = b-a
    w = [p+t for t in gaps]
    n = prod(w)
    r = n-128
    for q in w+[r]:
        certify_prime(q)
    assert p > n+1-r and 2*r > n+1
    assert sum(Fraction(a, t-a) for t in gaps[2:]) < 1
    assert sum(Fraction(b+s, t-b) for t in gaps[3:]) < 1
    e = len(w)-1
    uniform = p//(10*e*3**e*b*b*prod(gaps[3:]))
    assert uniform == 102
    pairs = list(combinations(w[:3], 2))
    first_miss = 1
    while all(pair_test(n*first_miss, x, y) for x, y in pairs):
        first_miss += 1
    assert first_miss == 12564
    print(f"Five factors, n = {n}, fourth certificate prime = {r}.")
    print(f"Uniform threshold through d={uniform}, exact pair test through d={first_miss-1}.")
    print(f"At d={first_miss} the pairs missing the test are "
          f"{[(x, y) for x, y in pairs if not pair_test(n*first_miss, x, y)]}.")
    for d in range(1, 120664):
        assert all(part_capacity(d*n, x, y) == 1 for x, y in pairs), d
    assert part_capacity(120664*n, w[0], w[2]) is None
    print("The position bounds certify all three pairs through d=120663.")
    return n, p, 120663


def minimum_for_pair(job):
    product, p, d = job
    n = d*product
    r = n-1 if n % 2 == 0 else n-2
    while not is_prime(r):
        r -= 2
    assert 2*r > n+1 and n+1-r <= p
    sizes = []
    for weight in (n, n+1):
        universal = universal_prime(weight)
        # Independent integer-bisection enumeration includes composite bases.
        witnesses = {q for q in candidate_bases(weight) if is_prime(q)}
        assert bool(universal) == bool(witnesses)
        if universal:
            assert universal in witnesses and sum(digits(weight, universal)) <= 2
        sizes.append(1 if universal else 4)
    return tuple(sizes), (n-r, d)


if __name__ == "__main__":
    small_controls()
    maximum_parts_controls()
    quartic_degree()
    product, p, last = five_factors()
    counts = Counter()
    largest_gap = (0, 0)
    with ProcessPoolExecutor(max_workers=4) as workers:
        jobs = ((product, p, d) for d in range(1, last+1))
        for sizes, gap in workers.map(minimum_for_pair, jobs, chunksize=16):
            counts[sizes] += 1
            largest_gap = max(largest_gap, gap)
    assert sum(counts.values()) == last
    print(f"Minimum sizes: {dict(sorted(counts.items()))}. Largest (gap, d): {largest_gap}.")
