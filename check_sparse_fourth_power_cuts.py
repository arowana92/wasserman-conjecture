"""Ternary digit-sum bounds and rational-polynomial digit gaps at cofactor 35."""

from math import comb, gcd


def first_carry(scale):
    """One-based position of the first carry in 16w + 19w, zero if none."""
    x, y = 16 * scale, 19 * scale
    position = 0
    while x and y:
        position += 1
        x, low_x = divmod(x, 3)
        y, low_y = divmod(y, 3)
        if low_x + low_y >= 3:
            return position
    return 0


def valuation(scale):
    x, y = 16 * scale, 19 * scale
    total = carry = 0
    while x or y:
        x, a = divmod(x, 3)
        y, b = divmod(y, 3)
        carry = (a + b + carry) // 3
        total += carry
    return total


def sparse_units(function, budget, stop_at=None):
    """Search every unit prefix with the given input digit-sum budget."""
    stop_at = {} if stop_at is None else stop_at
    layer = [(0, 0)]
    modulus = 1
    length = visited = max_failure = 0
    good = set()
    while layer:
        children = []
        visited += len(layer)
        for prefix, weight in layer:
            if weight == budget:
                failure = first_carry(function(prefix))
                max_failure = max(max_failure, failure)
                if not failure:
                    good.add(prefix)
                continue
            if prefix and 35 * function(prefix) < modulus:
                # All remaining zero digits work. Any pruning here needs
                # a separately proved obstruction to a further nonzero digit.
                assert prefix in stop_at, (prefix, weight, length)
                good.add(prefix)
                if length >= stop_at[prefix]:
                    continue
            for digit in range(int(length == 0), min(2, budget - weight) + 1):
                v = prefix + digit * modulus
                w = function(v)
                m = 3 * modulus
                if 16*w % m + 19*w % m < m:
                    children.append((v, weight + digit))
        length += 1
        modulus *= 3
        layer = children
    return sorted(good), length, visited, max(length, max_failure)


def offset_check(limit):
    stationary = []
    short = []
    for d in range(-limit, limit + 1):
        if not d % 2 or not d % 3:
            continue
        bound = 1 + max(35 * comb(4, i) * abs(d)**(4-i) for i in range(5))
        n, block = 1, 3
        while block <= bound:
            n += 1
            block *= 3
        for j in range(1, n + 4):
            s = (3**j + d) // 2
            if s <= 0 or not s % 2:
                continue
            v = valuation(s**4)
            if j < n:
                if v == 0:
                    short.append((d, j, s))
                continue
            next_v = valuation(((3**(j+4) + d)//2)**4)
            assert next_v >= v
            if v == 0:
                assert next_v == 0
                stationary.append((d, j % 4))
    assert stationary == [(1, 0)]
    assert {s for d, j, s in short} == {3281}
    assert len(short) == 7
    return stationary, short


def small_prefix_audit():
    # Independent low-to-high search, with arbitrary jumps between nonzero
    # digits. A carry of the current integer at layer p^k forces the next
    # nonzero digit to lie strictly below that modulus.
    good = set()
    stack = [(1, 1, 3), (2, 2, 3)]
    while stack:
        s, weight, place = stack.pop()
        w = s**4
        power = 1
        while power <= 35*w:
            power *= 3
            if 16*w // power + 19*w // power != 35*w // power:
                break
        else:
            good.add(s)
            if weight == 13:
                continue
            assert s == 3281
            power = 3**40
        while place < power and weight < 13:
            for digit in range(1, min(2, 13 - weight) + 1):
                root = s + digit * place
                # All lower layers were clear; check the changed column.
                v = root**4
                m = 3 * place
                if 16*v // m + 19*v // m == 35*v // m:
                    stack.append((root, weight + digit, m))
            place *= 3
    return sorted(good)


def high_digit_bound():
    # Test all leading prefixes with at most eight units. An arbitrary
    # fractional continuation is allowed at the right end of each prefix.
    layer = [(1, 1), (2, 2)]
    length = 1
    visited = 0
    while layer:
        visited += len(layer)
        children = []
        for h, weight in layer:
            a, b = h**4, (h+1)**4
            scale = 3**(3*length-3)
            while scale < 19*(b-a):
                scale *= 3
            possible = False
            for x in range(16*a//scale, (16*b-1)//scale + 1):
                for y in range(19*a//scale, (19*b-1)//scale + 1):
                    left, right = x, y
                    while left and right and left%3 + right%3 < 3:
                        left //= 3
                        right //= 3
                    if not left or not right:
                        possible = True
            if possible:
                for d in range(min(2, 8-weight) + 1):
                    children.append((3*h+d, weight+d))
        layer = children
        length += 1
    return length-1, visited


def bounded_sparse_search(budget, length_limit):
    # All successful roots within the digit and length bounds, including ones
    # shorter than the bound, are terminal inputs or zero-padded inputs.
    roots = set()
    visited = 0
    stack = [(1, 1, 3, 1), (2, 2, 3, 1)]
    while stack:
        s, weight, place, length = stack.pop()
        visited += 1
        if weight == budget or length == length_limit:
            if not first_carry(s**4):
                roots.add(s)
            continue
        m = 3*place
        for d in range(min(2, budget-weight)+1):
            t = s+d*place
            w = t**4
            if 16*w % m + 19*w % m < m:
                stack.append((t, weight+d, m, length+1))
    return sorted(roots), visited


def gap_checks():
    count = 0
    for p in (2, 3, 5):
        for denominator in (1, 2, 5, 7, 8):
            if gcd(p, denominator) != 1:
                continue
            period = 1
            while pow(p, period, denominator) != 1 % denominator:
                period += 1
            # The coefficient bound is below 64, including all row totals.
            # Eleven blanks suffice in every tested base and denominator.
            for g in range(12, 15):
                for h in (1, 12, 13, 14):
                    q, r = p**g, p**(g+h)
                    rows = []
                    for a, b, c in ((-7, -11, 7), (4, 6, 3), (-9, 12, 8)):
                        shift = -(a + b*q + c*r) % denominator
                        rows.append((a + shift, b, c))

                    def exact(left, right):
                        xs = [a + b*p**left + c*p**(left+right) for a, b, c in rows]
                        assert all(x > 0 and x % denominator == 0 for x in xs)
                        parts = [x // denominator for x in xs]
                        total = sum(parts)
                        power, out = p, 0
                        while power <= total:
                            out += total//power - sum(x//power for x in parts)
                            power *= p
                        return out

                    start = exact(g, h)
                    dg = exact(g + period, h) - start
                    assert dg >= 0 and exact(g + 3*period, h) == start + 3*dg
                    if h >= 12:
                        dh = exact(g, h + period) - start
                        assert dh >= 0 and exact(g + 2*period, h + 3*period) == start + 2*dg + 3*dh
                    count += 1
    return count


if __name__ == "__main__":
    print("Signed rational gap checks:", gap_checks())
    leaf = 3281
    derivative = sparse_units(lambda u: 4 * leaf**3 * u, 4)
    assert derivative == ([92], 6, 8, 6)
    assert min(3**40 - 35 * comb(4, j) * leaf**(4-j) * 92**j for j in range(5)) > 0
    assert first_carry(92**4) == 2
    table = sparse_units(lambda s: s**4, 13, {leaf: 40})
    assert table == ([3281, 265721], 70, 389454, 80)
    assert small_prefix_audit() == table[0]
    assert 3281 == 17*193 and 265721 == 41*6481
    assert gcd(17,193) == gcd(41,6481) == 1
    assert sparse_units(lambda s: s**4, 8) == ([], 39, 2620, 39)
    assert high_digit_bound() == (42, 14385)
    extended, visits = bounded_sparse_search(17, 80)
    assert extended == [3281, 265721, 21523361] and visits == 21703037
    print("Seventeen-unit roots and node count:", extended, visits)
    stationary, short = offset_check(10000)
    print("Derivative and fourth-power searches:", derivative, table)
    print("Stationary (offset, exponent mod 4):", stationary)
    print("Small exceptional representations:", short)
