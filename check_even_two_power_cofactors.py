"""Unbounded prime-power partners for cofactors D equal to 2^j times an odd prime power."""
from argparse import ArgumentParser
from collections import Counter, deque
from math import gcd, isqrt

from carry_blocks import CarryGraph, accepted_hull


def primes_and_powers(limit):
    primes = []
    composite = bytearray(limit+1)
    for p in range(3,limit+1,2):
        if composite[p]:
            continue
        primes.append(p)
        if p*p <= limit:
            composite[p*p::2*p] = bytes([1])*((limit-p*p)//(2*p)+1)
    powers = []
    for p in primes:
        power = p
        while power <= limit:
            powers.append((power,p))
            power *= p
    return primes, sorted(powers)


def merge_progressions(x,m,y,n):
    if m == 0:
        return (x,0) if (x-y)%n == 0 else None
    g = gcd(m,n)
    if (y-x)%g:
        return None
    period = m*(n//g)
    residue = x+m*((y-x)//g*pow(m//g,-1,n//g)%(n//g))
    return residue%period, period


def digit_subset(x,total,base):
    while x:
        if x%base > total%base:
            return False
        x //= base
        total //= base
    return True


def positions(w, powers, primes, limit):
    for P,p in powers:
        if w*P > limit:
            break
        for u in primes:
            if u >= min(w,P):
                break
            if P%u == 0:
                continue
            for a in range(u,w,2*u):
                if a == w-1 == u:
                    continue
                for b in range(u,((w-a)*P-1)//w+1,u):
                    if b%p == 0:
                        continue
                    r = pow(a*P,-1,w)
                    r += w*((-pow(w*P,-1,u)-r)*pow(w,-1,u)%u)
                    r += w*u*((pow(b*w,-1,P)-r)*pow(w*u,-1,P)%P)
                    slopes = a*P,b*w,(w-a)*P-b*w
                    yield P,p,u,a,b,r,w*u*P,slopes


def prime_power_possible(fixed, step, slopes):
    if step == 0:
        if fixed < 2:
            return False
        common = fixed
    else:
        common = gcd(fixed,step)
        if common <= 1:
            return True
    # A found divisor is enough to identify the only possible prime base.
    # If the trial budget does not factor this number, leave the case open.
    q = next((i for i in range(2,min(isqrt(common),10000)+1) if common%i == 0), None)
    if q is None:
        if common > 10000**2:
            return True
        q = common
    while common%q == 0:
        common //= q
    if common != 1:
        return False
    entries = slopes
    while any(entries):
        if sum(a%q for a in entries) >= q:
            return False
        entries = tuple(a//q for a in entries)
    return True


def prefix_graph(base, slopes, power, residue):
    low = residue%power
    values = [a*low+b for a,b in zip(slopes,(0,0,1))]
    digits = [x%power for x in values]
    while any(digits):
        if sum(x%base for x in digits) >= base:
            return None
        digits = [x//base for x in digits]
    return CarryGraph(base, slopes, [x//power for x in values])


def possible_suffix(graph, modulus, target):
    # For moderate moduli, enumerate residues of all accepted suffixes. Before
    # inserting a very large number of pairs, conservatively give up this test.
    parents = [[] for _ in graph.edges]
    for i,row in enumerate(graph.edges):
        for digit,j in row:
            parents[j].append((i,digit))
    reached = {(i,0) for i in range(len(parents)) if graph.accept >> i & 1}
    todo = list(reached)
    while todo:
        i,s = todo.pop()
        if i == 0 and s == target:
            return True
        for j,digit in parents[i]:
            key = j,(graph.base*s+digit)%modulus
            if key not in reached:
                reached.add(key)
                todo.append(key)
        if len(reached) > 60000:
            return True
    return False


def direct_progression(base, slopes, residue, modulus):
    # Read the digits of k in Q = residue + modulus*k. The direct traversal
    # can quickly exhaust graphs even when the suffix-residue product is large.
    initial = tuple(a*residue+b for a,b in zip(slopes,(0,0,1)))
    jumps = tuple(a*modulus for a in slopes)
    todo = deque([initial])
    seen = {initial}
    budget = max(100,120000//base)
    while todo:
        state = todo.popleft()
        terminal = state
        while any(terminal) and sum(x%base for x in terminal) < base:
            terminal = tuple(x//base for x in terminal)
        if not any(terminal):
            return True
        for digit in range(base):
            vals = tuple(x+a*digit for x,a in zip(state,jumps))
            if sum(x%base for x in vals) >= base:
                continue
            child = tuple(x//base for x in vals)
            if child not in seen:
                if len(seen) >= budget:
                    return True
                seen.add(child)
                todo.append(child)
    return False


def verify_shape(w,P,p,u,a,b,r,modulus,slopes):
    if not digit_subset(a*r%P,w*r%P,p):
        return 'fixed digit'
    if (a//u*P*r)%u + (b//u*w*r)%u >= u:
        return 'fixed digit'
    graph = prefix_graph(2,slopes,w,r)
    if graph is None:
        return 'binary prefix'
    hull = accepted_hull(graph)
    if hull is None:
        return 'empty binary language'
    h,d = hull
    combined = merge_progressions(r%w+w*h,w*d,r,modulus)
    if combined is None:
        return 'binary congruence'
    if not prime_power_possible(*combined,slopes):
        return 'prime-power obstruction'
    for base, power in ((2,w), (u,u), (p,P)):
        if not direct_progression(base,slopes,r,modulus):
            return 'progression graph'
        # High prime bases can have dense graphs. Skipping this optional test
        # only retains an extra shape and never certifies a weight.
        if base > 127:
            continue
        if base != 2:
            graph = prefix_graph(base,slopes,power,r)
        if graph is None or not possible_suffix(graph,modulus//power,r//power):
            return 'suffix residues'
    return None


def enumerate_cofactors(bound, small_multiplier_bound):
    limit = max(bound, small_multiplier_bound)
    primes, powers = primes_and_powers(limit//8)
    total = Counter()
    survivors = []
    w = 8
    while w*5 <= limit:
        local_bound = max(bound, small_multiplier_bound) if w in (8,16) else bound
        for P,p,u,a,b,r,modulus,slopes in positions(w,powers,primes,local_bound):
            reason = verify_shape(w,P,p,u,a,b,r,modulus,slopes)
            total[reason] += 1
            if reason is None:
                survivors.append((w,P,u,a,b))
                print("Unresolved shape",survivors[-1],flush=True)
        w *= 2
    print("Tested",sum(total.values()),"assignments.", dict(total),flush=True)
    return survivors


if __name__ == '__main__':
    cli = ArgumentParser(description=__doc__)
    cli.add_argument('--bound', type=int, default=15423)
    cli.add_argument('--small-multiplier-bound', type=int, default=40000)
    args = cli.parse_args()
    assert not enumerate_cofactors(args.bound,args.small_multiplier_bound)
