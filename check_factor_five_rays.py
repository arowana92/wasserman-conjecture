"""Check every binary exponent for the factor-five rows at bounded odd prime powers."""

from argparse import ArgumentParser
from collections import Counter, deque
from math import gcd, isqrt

from search_even_endpoint import prime_powers


def binary_prefixes(bits):
    def extend(index, value, previous, carry):
        if index == bits:
            if not (previous and (previous + carry) % 2):
                yield value
            return
        for digit in (0, 1):
            total = digit + previous + carry
            if not (previous and total % 2):
                yield from extend(index+1, value+(digit << index), digit, total//2)
    yield from extend(1, 1, 1, 0)


def exponent_classes(q, prime, first, counts):
    order = q - q//prime
    for factor in prime_powers(order):
        while order % factor == 0 and pow(16, order//factor, q) == 1:
            order //= factor

    step = isqrt(order)+1
    baby = {pow(16, j, q): j for j in range(step)}
    giant = pow(16, -step, q)
    z = 2*q
    while z:
        z = (z-1) & (2*q)
        if not z or z % 15 or gcd(z, q) > 1:
            continue
        counts['binary subsets'] += 1
        w = pow(z, -1, q)
        power = prime
        while power <= q and (3*w) % power <= (5*w) % power:
            power *= prime
        if power <= q:
            continue
        target = w*pow(2, -first, q) % q
        if pow(target, order, q) != 1:
            continue
        i = 0
        while target not in baby:
            target = target*giant % q
            i += 1
            if i > step:
                raise ArithmeticError('subgroup membership or discrete logarithm')
        k = first + 4*((i*step+baby[target]) % order)
        assert pow(2, k, q) == w
        yield z, k, 4*order


def no_carries(values, p):
    a,b,c = values
    while a or b or c:
        if a%p + b%p + c%p >= p:
            return False
        a,b,c = a//p,b//p,c//p
    return True


def initial_states(q, z, p, k, stride):
    modulus = p
    while pow(2, stride, modulus*p) == 1:
        modulus *= p
    if p == 3:
        modulus = max(modulus, 9)
    first = pow(2, k, modulus)
    ratio = pow(2, stride, modulus)
    w = first
    while True:
        parts = (3*q*w, z*w, (2*q-z)*w+1)
        if (p != 3 or (5*q*w+1) % 9) and no_carries([a%modulus for a in parts],p):
            yield tuple(a//modulus for a in parts)
        w = w*ratio % modulus
        if w == first:
            break


def successors(state, slopes, p):
    for digit in range(p):
        values = tuple(s+a*digit for s,a in zip(state,slopes))
        if sum(x%p for x in values) < p:
            yield tuple(x//p for x in values)


def accepts_scalar(q, z, p, k, stride):
    slopes = 3*q, z, 2*q-z
    pending = deque(initial_states(q,z,p,k,stride))
    visited = set(pending)
    while pending:
        state = pending.popleft()
        if no_carries(state,p):
            return True
        for successor in successors(state,slopes,p):
            if successor not in visited:
                visited.add(successor)
                pending.append(successor)
    return False


def search(bits):
    counts = Counter()
    survivors = []
    for q in binary_prefixes(bits):
        if q <= 5 or q%8 != 3 or bin(q).count('1') < 7:
            continue
        first = next((k for k in range(4) if q*(1<<k)%15 == 7), None)
        if first is None:
            continue
        factors = prime_powers(q)
        if len(factors) != 1:
            continue
        prime = next(iter(factors))
        if prime <= 5:
            continue
        counts['prime powers'] += 1
        for z,k,stride in exponent_classes(q,prime,first,counts):
            counts['classes'] += 1
            for p in (3,5):
                if not accepts_scalar(q,z,p,k,stride):
                    counts[f'excluded at {p}'] += 1
                    break
            else:
                survivors.append((q,z,k,stride))
    return counts, sorted(survivors)


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('bits', type=int, default=21, nargs='?')
    bits = parser.parse_args().bits
    if bits < 2:
        parser.error('supply at least two bits')
    counts, survivors = search(bits)
    print(dict(counts))
    for row in survivors:
        print('survivor', *row)
    first_ray = (1467779,2918400,2337019,2935556)
    second_ray = (2823523,5636100,3948614,5647044)
    if bits in (21,23):
        expected = (122,5840,1194,1173,20) if bits==21 else (269,20472,3912,3860,50)
        keys = ('prime powers','binary subsets','classes','excluded at 3','excluded at 5')
        assert counts == dict(zip(keys,expected))
        assert survivors == ([first_ray] if bits==21 else [first_ray,second_ray])
