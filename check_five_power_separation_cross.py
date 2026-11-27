"""Independent binary subset graph, rational reconstruction and finite range checks."""

from collections import defaultdict
from itertools import combinations
from math import gcd

from check_five_power_separation import SHAPES


def subset_graph(shape):
    k, c, q = shape
    B = 2 << k
    initial = (q, (c*q-1)//B, q >> k, (q+(2*q) % B)//B)
    states = [initial]
    index = {initial: 0}
    edges, reverse = [], [[]]
    stops = set()
    for s, (u, v, previous, carry) in enumerate(states):
        if v & ~u == 0 and not (previous and (previous+carry) % 2):
            stops.add(s)
        out = []
        for bit in (0,1):
            digit_sum = bit+previous+carry
            z = v+c*bit
            if (previous and digit_sum % 2) or (z % 2 > u % 2):
                continue
            next_state = (u//2+(bit << k), z//2, bit, digit_sum//2)
            if next_state not in index:
                index[next_state] = len(states)
                states.append(next_state)
                reverse.append([])
            t = index[next_state]
            reverse[t].append(s)
            out.append((bit,t))
        edges.append(out)
    good = set(stops)
    pending = list(stops)
    for t in pending:
        for s in reverse[t]:
            if s not in good:
                good.add(s)
                pending.append(s)
    return {s:[(b,t) for b,t in edges[s] if t in good] for s in good}, stops


def possible_low_words(shape, edges, bits):
    k, _, q = shape
    front = {q:0}
    for pos in range(k+1,bits):
        front = {v+(digit << pos):t for v,s in front.items() for digit,t in edges[s]}
    return set(front)


def logarithms(values, bits):
    # 11 also generates all classes 1 or 3 mod 8. Divide out successive powers.
    modulus = 1 << bits
    inverse = pow(11,-1,modulus)
    steps = [pow(inverse, 1 << j, modulus) for j in range(bits-2)]
    exponents = {}
    for v in values:
        word, exponent = v*inverse % modulus, 1
        for bit in range(1,bits-2):
            if word & (1 << (bit+2)):
                word = word*steps[bit] % modulus
                exponent += 1 << bit
        assert word == 1 and pow(11,exponent,modulus) == v
        exponents[v] = exponent
    return exponents


def short_fraction(ratio, modulus, bound):
    # Positive residuals from the continued-fraction convergents of ratio/modulus.
    previous, remainder = modulus, ratio
    prev_den, den = 0, 1
    while remainder and abs(den) <= bound:
        if den > 0 and remainder <= bound and den*remainder % 2 and gcd(den,remainder) == 1:
            return den, remainder
        quotient, next_rem = divmod(previous, remainder)
        previous, remainder = remainder, next_rem
        prev_den, den = den, prev_den-quotient*den
    return None


def check_separation(all_words):
    small = 1 << 124
    unit_order = small//4
    whole = 1 << 176
    by_residue = defaultdict(list)
    for words in all_words:
        for value in words:
            by_residue[value % small].append(value)
    logs = logarithms(by_residue, 124)
    inverse = {x:pow(i,-1,unit_order) for x,i in logs.items()}
    possible = 0
    for x,y in combinations(by_residue, 2):
        f = short_fraction(logs[y]*inverse[x] % unit_order, unit_order, 10**18)
        if f is None:
            continue
        a,b = f
        possible += 1
        # Power the actual binary prefixes directly at 176 bits.
        residues = {pow(z,b,whole) for z in by_residue[x]}
        assert all(pow(w,a,whole) not in residues for w in by_residue[y])
    assert possible == 3947
    print("Continued fractions recovered 3947 coarse relations. Direct powers excluded all.")


def digital_test(n, a, b, p):
    while n or a or b:
        if a % p + b % p > n % p:
            return False
        n, a, b = n//p, a//p, b//p
    return True


def check_scalar(shape, edges, stops, bits, expected):
    # Traverse only edges that allow precisely the remaining length. Compress
    # forced segments between branch points before writing out their digits.
    k, c, q = shape
    P = 1 << k
    length = bits-k-1
    can = [stops]
    for j in range(1,length+1):
        can.append({s for s,out in edges.items() if any(t in can[j-1] for _,t in out)})

    from functools import lru_cache
    @lru_cache(maxsize=None)
    def segment(s, remaining):
        suffix = offset = 0
        while remaining:
            choices = [(b,t) for b,t in edges[s] if t in can[remaining-1]]
            if len(choices) != 1:
                return remaining, offset, suffix, choices
            b,s = choices[0]
            suffix += b << offset
            offset += 1
            remaining -= 1
        return 0, offset, suffix, []

    # A different low-layer sieve from the first verifier.
    mod = 27*625
    allowed = bytearray(mod)
    for v in range(mod):
        if P*v % 15 != 7 or (5*P*v+1) % 9 == 0:
            continue
        n, a, b = 1+5*P*v, 3*P*v, c*v
        allowed[v] = all(digital_test(n % m, a % m, b % m, p) for p,m in [(3,27),(5,625)])
    count = 0
    todo = [(0,length,q)]
    while todo:
        s, remaining, value = todo.pop()
        j = bits-remaining
        rest, offset, suffix, choices = segment(s,remaining)
        value += suffix << j
        if rest == 0:
            count += 1
            if allowed[value % mod]:
                n, a, b = 1+5*P*value, 3*P*value, c*value
                assert not (digital_test(n,a,b,3) and digital_test(n,a,b,5))
        else:
            for b,t in choices:
                todo.append((t,rest-1,value+(b << (j+offset))))
    assert count == expected
    print("Trimmed binary traversal:", k, bits, count)


def check_primes(all_words):
    from math import isqrt
    def primes(bound):
        smallest = bytearray(bound+1)
        for d in range(2,isqrt(bound)+1):
            if not smallest[d]:
                for m in range(d*d,bound+1,d):
                    smallest[m] = 1
        return [p for p in range(11,bound+1,8) if not smallest[p]]

    mod = 1 << 130
    period = mod//4
    prime_logs = logarithms(primes(1_000_000), 130)
    inverse = {p:pow(x,-1,period) for p,x in prime_logs.items()}
    for shape, words, expected in zip(SHAPES,all_words,(8577,1275530,27148)):
        k,c,_ = shape
        P = 1 << k
        entries = logarithms({x % mod for x in words},130).values()
        count = 0
        for p,i in inverse.items():
            if not digital_test(5*P,3*P,c,p):
                continue
            for x in entries:
                e = x*i % period
                count += 1
                q = pow(p,e % 12,45)
                if P*q % 15 != 7 or (1+5*P*q) % 9 == 0:
                    continue
                for base in (5,3):
                    precision = base**32
                    residue = pow(p,e,precision)
                    n, a, b = (1+5*P*residue) % precision, (3*P*residue) % precision, (c*residue) % precision
                    if not digital_test(n,a,b,base):
                        break
                else:
                    raise AssertionError((k,p,e))
        assert count == expected
        print("Independent prime-power exponent count:", k, count)


if __name__ == "__main__":
    graphs = [subset_graph(shape) for shape in SHAPES]
    words = [possible_low_words(shape,edges,176) for shape,(edges,_) in zip(SHAPES,graphs)]
    check_separation(words)
    check_primes(words)
    for shape,(edges,stops),bits,count in zip(SHAPES,graphs,(1152,384,384),(24845,5321850,430247)):
        check_scalar(shape,edges,stops,bits,count)
