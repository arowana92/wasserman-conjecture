"""Audit the binary shapes at two upper powers and count their possible offsets."""

from collections import Counter

from check_dense_offsets import free_rows, has_size_obstruction


def pattern_counts():
    block_counts = {}
    for h in range(2, 12):
        mask = (1 << (h+1))-1
        forbidden = 1 | (1 << (h-1)) | (1 << h)
        counts = {0: 1}
        for _ in range(128):
            successor = Counter()
            for tail, count in counts.items():
                shifted = (tail << 1) & mask
                successor[shifted] += count
                if not tail & forbidden:
                    successor[shifted+1] += count
            counts = successor
        block_counts[h] = sum(counts.values())
    masks = [x for x in range(1 << 12) if not (x & (x << 1))]
    grids = {x: 1 for x in masks}
    for _ in range(4):
        grids = {upper: sum(count for lower, count in grids.items()
                            if not upper & (lower | (lower << 1))) for upper in masks}
    grid_count = sum(grids.values())
    assert max(block_counts.values()) < 2**64
    assert grid_count == 3396136349 and grid_count**100 < 2**3180
    print('Largest 128-bit block count:', max(block_counts.values()))
    print('Five by twelve grid patterns:', grid_count)


def binary_support_check(limit):
    # Enumerate every scalar directly and add the shifted copies as bit masks.
    direct = set()
    for h in range(1, limit.bit_length()):
        r = 1 << h
        for u in range(1, limit//(3*r+1)+1):
            d = (3*r+1)*u
            if (u | r*u | 2*r*u) == d:
                direct.add((d, h))

    # Independently build bit sets one position at a time with forbidden distances.
    recovered = set()
    for h in range(1, limit.bit_length()):
        prefix = (3 << h)+1
        maximum = limit//prefix
        words = [0]
        for bit in range(maximum.bit_length()):
            unit = 1 << bit
            blockers = (unit >> 1) | (unit >> h) | (unit >> (h+1))
            words += [word+unit for word in words if not (word & blockers) and word+unit <= maximum]
        recovered.update((prefix*word, h) for word in words if word)
    assert recovered == direct
    offsets = {d for d, _ in direct}
    assert all(bin(d).count('1')%3 == 0 for d in offsets)
    print('Candidate offsets through', limit, ':', len(offsets), 'distinct;', len(direct), 'offset and gap pairs')
    print('First candidates:', sorted(offsets)[:16])
    print('Candidates with gap at least two:', len({d for d, h in direct if h >= 2}))


def check_actual_rows():
    seen = survivors = cases = 0
    for k in (79, 103, 151, 203):
        q = 1 << k
        for h in (1, 2, 3, 4, 5, 6, 7, 8, k//2, k-1):
            r, middle = 1 << h, 1 << (k-h)
            for d in range(1, 257):
                if d**25 > q:
                    continue
                n = 3*q+middle+d
                cases += 1
                for row in free_rows(n, 2):
                    seen += 1
                    if has_size_obstruction(n, row):
                        continue
                    u, rem = divmod(d, 3*r+1)
                    assert rem == 0 and (u | r*u | 2*r*u) == d
                    f = middle+u
                    assert row == (f, r*f, 2*r*f)
                    survivors += 1
    assert survivors
    print('Weights, binary-free rows, rows passing the size test:', cases, seen, survivors)


if __name__ == '__main__':
    pattern_counts()
    binary_support_check(1000000)
    check_actual_rows()
