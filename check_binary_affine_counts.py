"""Test prefix bounds against known digits, factorial floors and sharp binary examples."""

from itertools import combinations_with_replacement
from math import comb

from carry_blocks import CarryGraph
from check_three_factors import valuation


def order(value, base):
    exponent = 0
    while value % base == 0:
        exponent += 1
        value //= base
    return exponent


def prefix_test(base, slopes, offsets, length):
    # Each product for t mod p^L is known modulo p^(L+v_p(slope)).
    slopes = (*slopes, sum(slopes))
    offsets = (*offsets, sum(offsets))
    known = [length + order(x, base) for x in slopes]
    places = [base**i for i in range(max(known))]
    columns = [[i for i, last in enumerate(known) if j < last]
               for j in range(max(known))]
    count = accepted = 0
    for t in range(base**length):
        row = [a*t + b for a, b in zip(slopes, offsets)]
        live = True
        for place, indices in zip(places, columns):
            used = 0
            for i in indices:
                digit = row[i] // place % base
                used += base-1-digit if i == len(row)-1 else digit
            if used > base-1 or (len(indices) == len(row) and used != base-1):
                live = False
                break
        count += live
        if valuation(row[-1], row[:-1], base) == 0:
            assert live, (t, slopes, offsets, length)
            accepted += 1
    delay = max(known)-min(known)
    if base == 2:
        # Certified upper rational approximations, using lambda^k-lambda^(k-1)=1.
        n, d = {2: (16181, 10000), 3: (14656, 10000), 4: (13803, 10000)}[len(row)-1]
        k = len(row)-1
        assert n**k-n**(k-1)*d > d**k
        assert count*d**(length+delay) <= n**(length+delay)
    else:
        assert count * 3**length <= 2**(delay*(base-1)) * (2*base+1)**length
    return count, accepted


def sharp_counts():
    for k in range(2, 7):
        graph = CarryGraph(2, tuple(2**i for i in range(k)), (0,)*k)
        # Count paths to states whose remaining multiplication carries can terminate.
        counts = [(graph.accept >> i) & 1 for i in range(len(graph.edges))]
        for length in range(129):
            # Compress the k-1 required zeros between successive ones.
            expected = 1 + sum(comb(length-(k-1)*(ones-1), ones)
                               for ones in range(1, (length+k-1)//k+1))
            assert counts[0] == expected, (k, length, counts[0], expected)
            counts = [sum(counts[j] for _, j in edges) for edges in graph.edges]


def main():
    sharp_counts()
    cases = prefixes = exact = 0
    for base, k, largest, length in [(2, 2, 9, 8), (2, 3, 8, 8), (2, 4, 6, 7),
                                     (3, 3, 7, 4), (5, 3, 7, 3), (7, 3, 6, 3)]:
        for slopes in combinations_with_replacement(range(1, largest+1), k):
            for offsets in [(0,)*k, (0,)*(k-1)+(1,), tuple(i*i+2 for i in range(k))]:
                n, m = prefix_test(base, slopes, offsets, length)
                prefixes += n
                exact += m
                cases += 1
    assert prefixes > exact > 0
    print(f'{cases} shape/base checks, {prefixes} admissible prefixes, '
          f'{exact} inputs confirmed by factorial floors.')
    print('Sharp counts matched for 2 to 6 summands at each bit length 0 through 128.')


if __name__ == '__main__':
    main()
