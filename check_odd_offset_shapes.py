"""Test the two necessary odd-prime shapes and count the eligible offsets they cover."""

from collections import Counter

from check_dense_offsets import free_rows, has_size_obstruction


def adds_freely(parts, p):
    numbers = list(parts)
    while any(numbers):
        digits, numbers = zip(*[(v%p, v//p) for v in numbers])
        if sum(digits) >= p:
            return False
    return True


def offset_lists(limit, p):
    equal, scaled = set(), set()
    powers = []
    power = p
    while power <= limit:
        powers.append(power)
        power *= p
    # Generate all possible values P+D in the digit alphabet 0,3,6,...
    totals = [0]
    place = 1
    while 3*place <= 2*limit:
        totals = [x+digit*place for x in totals for digit in range(0, p, 3)
                  if x+digit*place <= 2*limit]
        place *= p
    for total in totals:
        for middle in powers:
            d = total-middle
            if middle <= d <= limit and d%p == 0:
                equal.add(d)
    for r in powers:
        for u in range(1, limit//(3*r+1)+1):
            d = (3*r+1)*u
            if d%p not in (0, 1):
                continue
            if any(adds_freely((r*u, c*r*u, ((2-c)*r+1)*u), p) for c in (1, 2)):
                scaled.add(d)
    eligible = limit//p + (limit-1)//p + 1
    assert not any(d%(2*p) == 1 for d in equal | scaled)
    return eligible, equal, scaled


def direct_lists(limit, p):
    equal, scaled = set(), set()
    for d in range(1, limit+1):
        if d%p not in (0, 1):
            continue
        power = p
        while power <= d:
            total = d+power
            if d%p == 0 and total%3 == 0 and adds_freely((total//3,)*3, p):
                equal.add(d)
            u, remainder = divmod(d, 3*power+1)
            if remainder == 0:
                for c in (1, 2):
                    low = (power*u, c*power*u, ((2-c)*power+1)*u)
                    modulus = p
                    while modulus <= d and sum(x%modulus for x in low) == d%modulus:
                        modulus *= p
                    if modulus > d:
                        scaled.add(d)
            power *= p
    return equal, scaled


def sample_control(p):
    counts = Counter()
    for k in (55, 101):
        q = p**k
        for ell in sorted({1, 2, k//2, k-2, k-1}):
            middle = p**ell
            r = q//middle
            for d in sorted(set(range(1, 65)) | {3*p+1, 3*p*p+1}):
                if d*(d-1)%p or d**25 > q:
                    continue
                n = 3*q+middle+d
                counts['weights'] += 1
                for row in free_rows(n, p):
                    counts['rows'] += 1
                    if has_size_obstruction(n, row):
                        continue
                    if middle <= d:
                        assert row[0] == row[2] and d%p == 0
                        counts['equal survivors'] += 1
                    else:
                        u, rem = divmod(d, 3*r+1)
                        assert rem == 0 and u
                        f = middle+u
                        choices = [tuple(sorted((r*f, c*r*f, ((2-c)*r+1)*f))) for c in (1, 2)]
                        assert row in choices
                        counts['scaled survivors'] += 1
    assert counts['equal survivors'] and counts['scaled survivors']
    return counts


def check():
    for p in (5, 7, 11):
        for incoming in (0, 1):
            choices = [x for x in range(p) if x+(2*x+incoming)%p < p]
            assert len(choices) == (p+1)//2
        eligible, equal, scaled = offset_lists(1000000, p)
        assert ({d for d in equal if d <= 5000}, {d for d in scaled if d <= 5000}) == direct_lists(5000, p)
        print('Prime, eligible, equal, scaled, distinct candidates:', p, eligible,
              len(equal), len(scaled), len(equal | scaled))
        print('Direct sample:', dict(sample_control(p)))


if __name__ == '__main__':
    check()
