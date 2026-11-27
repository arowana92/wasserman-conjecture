"""Check the digit-gap, block orbit and quadratic-window obstructions."""
from math import comb, gcd, lcm


def digits(n,p):
    result = []
    while n:
        result.append(n % p)
        n //= p
    return result


def fits(r,t,p):
    while r:
        if r % p > t % p:
            return False
        r //= p
        t //= p
    return True


def blocks(t,p,d):
    occupied = [j for j,a in enumerate(digits(t,p)) if a]
    starts = [0]
    for left,right in zip(occupied, occupied[1:]):
        if right-left > d:
            starts.append(right)
    ends = starts[1:]+[occupied[-1]+1]
    return [(i,p**(j-i), t//p**i % p**(j-i)) for i,j in zip(starts, ends)]


def row_checks():
    hits = nontrivial = 0
    for p in (2,3,5,7):
        for n in range(1+p, 165, p):
            t = n-1
            while t % p == 0:
                t //= p
            w = (n-1)//t
            for a in range(1,n//3+1):
                for b in range(a,(n-a)//2+1):
                    c = n-a-b
                    if comb(n,a)*comb(n-a,b) % p == 0:
                        continue
                    part, = [v for v in (a,b,c) if v % w == 1]
                    s = t//gcd(part,t)
                    for j,digit in enumerate(digits(t,p)):
                        if j and digit:
                            power = p**j
                            assert s*(1+w*(t % power)) >= power
                            nontrivial += 1+w*(t % power) < power
                    hits += 1
    assert nontrivial
    print("Exact coefficient checks",hits,"layers requiring s > 1",nontrivial)


def orbit_checks():
    transitions = bounded_runs = sharp = clustered_pairs = 0
    for p in (2,3,5,7):
        for t in range(3, 651):
            if t % p == 0:
                continue
            orbit = [1]
            while (nxt := orbit[-1]*p % t) != 1:
                orbit.append(nxt)
            ok = [fits(r,t,p) for r in orbit]
            period = len(orbit)
            for step in (1,2,3):
                pieces = blocks(t,p,step)
                if len(pieces)<2:
                    continue
                M = min(max(lcm(a,b) for _,_,b in pieces) for _,_,a in pieces)
                limit = 0
                while p**(step*(limit+1)) <= M:
                    limit += 1
                for i in range(gcd(period,step)):
                    cycle = [(i+j*step) % period for j in range(period//gcd(period,step))]
                    assert not all(ok[j] for j in cycle)
                    run = best = 0
                    for j in cycle+cycle:
                        run = run+1 if ok[j] else 0
                        best = max(best,run)
                    assert best <= limit+1
                    bounded_runs += 1
                    sharp += best == limit+1 and best >= 2
                for j,r in enumerate(orbit):
                    n = (j+step) % period
                    if not (ok[j] and ok[n]):
                        continue
                    q = (p**step*r)//t
                    for i,modulus,a in pieces:
                        assert p**step*(r//p**i % modulus) == q*a+orbit[n]//p**i % modulus
                    transitions += 1
            for n in range(1,len(digits(t,p))):
                P = p**n
                H,L = divmod(t,P)
                bound = H*L
                det = [H*(r % P)-L*(r//P) for r in orbit]
                step = 1
                while (p**step+1)*bound < t:
                    for i,v in enumerate(det):
                        j = (i+step) % period
                        if ok[i] and ok[j]:
                            assert det[j] == p**step*v and v
                            clustered_pairs += 1
                    step += 1
    assert min(transitions,sharp,clustered_pairs)>0
    print("Block transitions",transitions,"orbit run bounds",bounded_runs,
          "attained run bounds of at least two",sharp,"short determinant gaps",clustered_pairs)


def initial_windows():
    checked = isolated = 0
    for p in (2,3,5,7):
        for h in (9,16,30,48):
            for H,L in ((1,3),(3,5),(4,7),(5,11),(11,7)):
                if gcd(p,H*L)>1:
                    continue
                P = p**h
                total = H*P+L
                for e in range(1,h*h+1):
                    if (2*H*L)**(2*(e+2*h)) >= P**h:
                        break
                    r = pow(pow(p,e,total),-1,total)
                    assert not fits(r,total,p)
                    checked += 1
                for e in range(1,h+1):
                    if 1+p**e*L < P:
                        r = pow(pow(p,e,total),-1,total)
                        assert not fits(r,total,p)
                        isolated += 1
    print("Quadratic-window inverses",checked,"individual gap inverses",isolated)


if __name__ == "__main__":
    row_checks()
    orbit_checks()
    initial_windows()
