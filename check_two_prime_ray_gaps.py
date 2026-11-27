"""Check all finite carry graphs used for the leading-digit bounds."""

from fractions import Fraction
from itertools import product

from check_binary_windows import ln_bounds


def carry_graph(slopes, p, box):
    states = list(product(*box))
    graph = {}
    for s in states:
        graph[s] = []
        for digit in range(p):
            v = [u+digit*a for u, a in zip(s, slopes)]
            if sum(w % p for w in v) < p:
                graph[s].append((digit, tuple(w//p for w in v)))
    end = (0, 0, 0)
    good = {end}
    while True:
        prev = {s for s in states if any(t in good for _, t in graph[s])}
        if prev <= good:
            break
        good.update(prev)
    depth = 0
    while True:
        after = {t for s in good for _, t in graph[s] if t in good}
        if after == good:
            break
        good = after
        depth += 1
    assert end in good
    return depth, {s: [(b, t) for b, t in graph[s] if t in good] for s in good}


def certificate(slopes, p, divisor, labels):
    # For a positive completion y at s, divisor*y + labels[s] is a power of p.
    box = [range(-1, a+1) for a in slopes]
    depth, graph = carry_graph(slopes, p, box)
    positive = set()
    while True:
        prev = {s for s in graph if any(b or t in positive for b, t in graph[s])}
        if prev == positive:
            break
        positive = prev
    assert positive == set(labels), (slopes, p, positive, labels)

    # A zero suffix can terminate only after feeding zero digits down to the origin.
    zero = {(0, 0, 0)}
    while True:
        prev = {s for s in graph if any(b == 0 and t in zero for b, t in graph[s])}
        if prev <= zero:
            break
        zero.update(prev)
    for s, offset in labels.items():
        assert 1 <= offset <= divisor
        for b, t in graph[s]:
            if t in zero and b:
                assert is_power(divisor*b+offset, p)
            if t in positive:
                assert divisor*b+offset == p*labels[t]
    print(slopes, p, 'discard', depth, 'core', len(graph), 'positive states', len(positive))
    return depth


def graphs():
    for p, d, labels, g in [
        (2, 3, {(0, 1, 2): 1, (0, 2, 4): 2}, 4),
        (3, 2, {(0, 1, 3): 1}, 3),
    ]:
        assert certificate((1, 3, 6), p, d, labels) == g
    for p, labels in [
        (2, {(0, 0, 1): 1, (1, 1, 3): 3, (1, 2, 4): 4, (0, 1, 2): 2}),
        (3, {(0, 0, 1): 1, (0, 1, 2): 2, (1, 2, 4): 4, (1, 1, 3): 3}),
    ]:
        assert certificate((2, 3, 5), p, 5, labels) == 3


def independent_suffixes():
    # Directly divide the complete affine parts and test digits, without graph
    # transitions. Include all bounded incoming carries and large sparse integers.
    for slopes in ((1, 3, 6), (2, 3, 5)):
        for p in (2, 3):
            d = (3 if p == 2 else 2) if slopes == (1, 3, 6) else 5
            depth = 4 if (p, slopes) == (2, (1, 3, 6)) else 3
            scale = p**depth
            inputs = set(range(1, 512))
            for j in (12, 25, 80, 201):
                inputs.update(p**j//d + i for i in range(-scale, scale+1))
            accepted = 0
            for v in inputs:
                lo = d*(v//scale)
                centers = [lo+i for i in range(1, d+1)]
                for c in product(*(range(-1, a+1) for a in slopes)):
                    parts = [a*v+b for a,b in zip(slopes,c)]
                    while any(parts) and sum(a % p for a in parts) < p:
                        parts = [a//p for a in parts]
                    if not any(parts):
                        assert any(is_power(n, p) for n in centers)
                        accepted += 1
            assert accepted
            print(slopes, p, 'directly accepted affine inputs', accepted)


def is_power(n, p):
    while n > 1 and n % p == 0:
        n //= p
    return n == 1


def logarithm_bounds():
    l2, u2 = ln_bounds(Fraction(2), 3500)
    l3, u3 = ln_bounds(Fraction(3), 3500)
    low, high = l2/u3, u2/l3
    for bound, expected in [(10**100, 57), (10**1000, 3310)]:
        a,b = low.as_integer_ratio()
        c,d = high.as_integer_ratio()
        prev, denom = 1, 0
        maximum = 1
        while denom <= bound:
            t, rem = divmod(a,b)
            u, other = divmod(c,d)
            assert t == u
            a,b,c,d = b,rem,d,other
            prev, denom = denom, prev+t*denom
            if prev:
                maximum = max(maximum, 1+(denom+prev-1)//prev)
        assert maximum == expected
        print(f'Continued-fraction spacing constant {maximum} through 10^{len(str(bound))-1}')


if __name__ == '__main__':
    graphs()
    independent_suffixes()
    logarithm_bounds()
