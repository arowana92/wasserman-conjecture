"""Certify forced divisors near (1,1,m) for every m+1 through a million."""
from collections import deque
from itertools import product
from math import comb, isqrt
from functools import lru_cache

LIMIT = 1_000_000


def is_prime(p):
    return p >= 2 and all(p % d for d in range(2, isqrt(p)+1))


def difference_of_powers(n, p):
    while n % p == 0:
        n //= p
    n += 1
    while n % p == 0:
        n //= p
    return n == 1


def all_small_blocks(p, limit):
    blocks = [0]
    place = 1
    while place <= limit:
        blocks += [x + d*place for d in range(1, min((p-1)//2, limit//place)+1)
                   for x in blocks if x+d*place <= limit]
        place *= p
    return [x for x in blocks if x % p]


def sparse_witness(m, p, blocks):
    for x in blocks:
        modulus = p-1
        while modulus < 2*x:
            modulus = p*(modulus+1)-1
        while modulus <= m*x:
            d, e = x, m*x % modulus
            while d and 2*(d%p) + e%p < p:
                d //= p
                e //= p
            if not d:
                return x, modulus
            modulus = p*(modulus+1)-1


def graph_witness(m, p):
    # Iterative depth-first traversal. A back edge outside the zero state
    # supplies the digits of a periodic multiplier, least significant first.
    assert (m+1) % p == 0
    factor = (m+1)//p
    size = m//2 + 1
    done = bytearray(size)
    done[0] = 1
    for start in range(1, size):
        if done[start]:
            continue
        path = [start]
        next_digit = [0]
        active = {start: 0}
        while path:
            i = path[-1]
            digit = next_digit[-1]
            if digit > min(i % p, p-1-i%p):
                done[i] = 1
                del active[i]
                path.pop()
                next_digit.pop()
                continue
            next_digit[-1] += 1
            j = i//p + factor*digit
            assert j < size
            if done[j]:
                continue
            if j in active:
                steps = [d-1 for d in next_digit[active[j]:]]
                x = sum(d*p**a for a,d in enumerate(steps))
                assert x
                return x, p**len(steps)-1
            active[j] = len(path)
            path.append(j)
            next_digit.append(0)


@lru_cache(maxsize=None)
def base_length(power, p):
    h = 0
    while power % p == 0:
        power //= p
        h += 1
    assert power == 1 and h > 0
    return h


def digit_sum(z, p):
    total = 0
    while z:
        z, digit = divmod(z,p)
        total += digit
    return total


def check_witness(m, p, x, modulus):
    # Independent verification by Legendre digit sums on two whole
    # repetitions. The block arithmetic makes the construction work at
    # every repetition count.
    base_length(modulus+1, p)
    j, y = divmod(m*x, modulus)
    assert 0 < x and 2*x+y <= modulus and j <= m//2
    repeat = modulus+2
    w = x*repeat
    third = m*w+j
    assert digit_sum(2*w+third,p) == 2*digit_sum(w,p)+digit_sum(third,p)


def exceptional_ranks():
    # Reverse search removes exactly the states whose every path reaches 0.
    m, p = 18197, 3
    size = m//2+1
    parents = [[] for _ in range(size)]
    degrees = [0]*size
    for s in range(1,size):
        for d in range(p):
            j, y = divmod(m*d+s,p)
            if 2*d+y < p:
                parents[j].append(s)
                degrees[s] += 1
    queue = deque([0])
    ranks = [0]*size
    while queue:
        j = queue.popleft()
        for s in parents[j]:
            ranks[s] = max(ranks[s], 1+ranks[j])
            degrees[s] -= 1
            if not degrees[s]:
                queue.append(s)
    assert not any(degrees) and max(ranks) == 72 and size == 9099
    for s in range(1,size):
        for d in range(1+min(s%3,2-s%3)):
            assert ranks[s] > ranks[s//3 + 6066*d]
    return max(ranks)


def census():
    positives = negatives = 0
    fallback = []
    for p in range(3,isqrt(LIMIT)+2,2):
        if not is_prime(p):
            continue
        blocks = all_small_blocks(p,80000)
        for n in range(p*(p-1),LIMIT+1,p*(p-1)):
            if difference_of_powers(n,p):
                negatives += 1
                continue
            proof = sparse_witness(n-1,p,blocks)
            if proof is None:
                proof = graph_witness(n-1,p)
                fallback.append((n,p,base_length(proof[1]+1,p) if proof else None))
            if proof:
                check_witness(n-1,p,*proof)
                positives += 1
            else:
                assert (n,p) == (18198,3)
    assert [(n,p) for n,p,h in fallback if h is None] == [(18198,3)]
    print('Nontrivial pairs with periodic witnesses:',positives)
    print('Odd-prime difference-of-powers pairs:',negatives)
    print('Full-graph cases (m+1, prime, witness length):',fallback)
    print('Longest path to zero at m=18197:',exceptional_ranks())


def leading_cut(slopes, p, valuation, denominator, depth):
    # Exact rational caps, retaining the carry from the unread suffix.
    total = sum(slopes)
    scale = total*denominator
    front = {(0,0,0,0,0)}
    for h in range(1,depth+1):
        following = set()
        for a,b,c,out,spent in front:
            if spent+out >= valuation:
                continue
            for i,j,k in product(range(p),repeat=3):
                for incoming in range(3):
                    digit = i+j+k+incoming-p*out
                    if not 0 <= digit < p or (h == 1 and not digit):
                        continue
                    parts = p*a+i, p*b+j, p*c+k
                    n = sum(parts)+incoming
                    if all(v*scale < (n+1)*(s*denominator+total) and
                           n*(s*denominator-total) < (v+1)*scale
                           for v,s in zip(parts,slopes)):
                        following.add((*parts,incoming,spent+out))
        front = following
    assert not front, (slopes,p,valuation,depth)


if __name__ == '__main__':
    census()
    for row in ((5,2,1,50,4),(5,3,1,200,4),
                (17,2,1,500,6),(17,3,2,10**7,11),
                (19,2,2,10000,10),(19,5,1,10**6,7),
                (23,2,3,10**6,14),(23,3,1,10**6,9)):
        m, *params = row
        leading_cut((1,1,m),*params)
    for row in (((1,2,4),3,1,100,3),((1,3,7),2,3,10**6,12),
                ((1,3,7),5,1,10000,5),((5,13,19),7,1,10000,4),
                ((5,13,19),2,1,10**6,13)):
        leading_cut(*row)
    for p in (2,3,5,7):
        for m in range(1,25):
            t = p**2
            n = (m+2)*t-1
            value = comb(n,t)*comb(n-t,t-1)
            a, b = value, m+1
            while b%p == 0:
                assert a%p == 0
                a //= p
                b //= p
            assert a%p
    print('All thirteen cone bounds and 96 exact sharpness controls passed.')
