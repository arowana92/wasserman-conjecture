"""Check the high-power theorem by bounding each symbolic allocation separately."""

from itertools import combinations, product
from math import gcd, lcm

from check_three_factors import powers, valuation


def binomial_cap(a, b, i, j, unit):
    # Return K,v with the relevant divisor <= K*R^v on H=R^j, J=R^i.
    determinant = abs((-a)**i-(-b)**j)
    if determinant:
        return determinant, 0
    assert i % 2 == j % 2 == 1 and j >= 3
    return (2 if unit else 2*a*b), 1


def projected_cap(c, i, j, tags, offsets):
    # Tag 1 takes H, tag 2 takes J, and tag 3 takes both.
    left, middle = tags[:2]
    a, b = offsets[:2]
    if left == 0 or middle == 0:
        common = ((a if left == 0 else b), 0)
    elif left == 1:
        common = binomial_cap(a, b, i, j, True)
    else:
        common = binomial_cap(b, a, i, j, True)
    k, v = common[0]**3, 3*common[1]
    for tag, s in zip((left, middle), (a, b)):
        if tag == 0:
            bound, degree = s, 0
        elif tag == 3:
            bound, degree = c-s, 0
        elif tag == 1:
            bound, degree = binomial_cap(s, c-s, i, j, False)
        else:
            bound, degree = binomial_cap(c-s, s, i, j, False)
        k *= bound**2
        v += 2*degree
    assert v < 2*j
    return k, 2*j-v


def coprime_part(n, exclude):
    while gcd(n, exclude) > 1:
        n //= gcd(n, exclude)
    return n


def check_offset(c):
    top = (2*c)**64
    tiles = []
    for first, second in product(range(3), repeat=2):
        tags = tuple(int(first == s)+2*int(second == s) for s in range(3))
        for a in range(c+1):
            for b in range(c-a+1):
                offsets = a, b, c-a-b
                if any(tag == s == 0 for tag, s in zip(tags, offsets)):
                    continue
                if any((tag == 0 and s < 3) or (tag in (1, 2) and s == 0)
                       for tag, s in zip(tags[:2], offsets[:2])):
                    continue
                tiles.append((tags, offsets))
    survivors = set()
    numeric = 0
    j = 2
    while (c+1)**j <= top:
        for i in range(1, j):
            if gcd(i, j) != 1:
                continue
            for tags, offsets in tiles:
                k, degree = projected_cap(c, i, j, tags, offsets)
                # Allow every integer R>c within the cap, including composites.
                r = c+1
                while 27*r**degree <= k and r**j <= top:
                    numeric += 1
                    high, low = r**j, r**i
                    units = 0, high, low, high+low
                    x, y, _ = [units[tag]+s for tag, s in zip(tags, offsets)]
                    d = high+low+c
                    scale = coprime_part(gcd(x, y), d*r)
                    if scale >= 3 and scale**3 * gcd(d, x*y)**2 >= 27*d*d:
                        # Include a full prime power of d only if it fits in one slot.
                        sx = coprime_part(d, d//gcd(d, x))
                        sy = coprime_part(d, d//gcd(d, y))
                        if scale**3 * lcm(sx, sy)**2 >= 27*d*d:
                            atoms = powers(r)
                            if len(atoms) == 1 and next(iter(atoms)) > c:
                                survivors.add((r, i, j, min(x,y), max(x,y), d, scale))
                    r += 1
        j += 1
    expected = {(13, 1, 2, 18, 171, 190, 9)} if c == 8 else set()
    assert survivors == expected, (c, survivors)
    return numeric


def small_bases(c):
    checked = 0
    for p in range(2, c+1):
        if powers(p) != {p: p}:
            continue
        pp = []
        v = p
        while v <= (2*c)**64:
            pp.append(v)
            v *= p
        # A constant partition can be checked by an exact factorial valuation.
        offsets = [(a,b,c-a-b) for a in range(c+1) for b in range(c-a+1)
                   if valuation(c, (a,b,c-a-b), p) == 0]
        for low, high in combinations(pp, 2):
            if high//low <= c:
                continue
            d = high+low+c
            if low <= c and 27*d*d >= (low+c)**5:
                continue
            checked += 1
            if low <= c:
                candidates = ((x, y) for x in range(1, d-1) for y in range(x, d-x)
                              if valuation(d, (x,y,d-x-y), p) == 0)
            else:
                candidates = ((a+low*(i == 0)+high*(j == 0), b+low*(i == 1)+high*(j == 1))
                              for a, b, _ in offsets for i, j in product(range(3), repeat=2))
            for x, y in candidates:
                if not 0 < x <= y or x+y >= d:
                    continue
                scale = gcd(x,y)
                if scale < 3 or scale**3*gcd(d,x*y)**2 < 27*d*d:
                    continue
                assert coprime_part(scale, p*d)**3*gcd(d,x*y)**2 < 27*d*d
    return checked


def rational_approximants():
    tests = resonances = zero_index = 0
    for n in range(2, 90):
        for h in range(1, n):
            _, i, j = min((abs(i*n-j*h), i, j) for j in range(1, 9) for i in range(j+1)
                          if gcd(i, j) == 1)
            assert 8*abs(i*n-j*h) <= n
            zero_index += i == 0
            if i == 0:
                continue
            # Check the divisibility before taking absolute values. This covers
            # both signs of the exponent difference, including a zero difference.
            for a, b in product(range(1, 5), repeat=2):
                d = gcd(3**n+a, 3**h+b)
                gap = i*n-j*h
                diff = (-b)**j * 3**max(0, gap) - (-a)**i * 3**max(0, -gap)
                assert diff % d == 0
                tests += 1
                resonances += diff == 0
    assert resonances and zero_index
    print(f'{tests} gcd congruences, {resonances} exact equalities, '
          f'{zero_index} small-exponent approximants.', flush=True)


if __name__ == '__main__':
    rational_approximants()
    for offset in range(3, 9):
        print(f'Offset {offset}: {check_offset(offset)} integer root and row checks, '
              f'{small_bases(offset)} small-base pairs.', flush=True)
