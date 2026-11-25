"""Check finite intersections of carry languages and certify their two-base counting bounds."""

from argparse import ArgumentParser
from math import comb

from carry_blocks import CarryGraph, forced_bases, intersection
from check_endpoint_pair_bound import prime_powers
from check_prime_power_cofactors import shapes, survives
from search_endpoint import atoms, smallest_factors

# Remaining odd-Q cofactors from the full finite-automaton sieve through 1000.
COFACTORS = (364, 380, 589, 629, 646, 703, 710, 780, 850, 874, 880, 946, 950, 962, 986)


def small_crosscheck():
    cases = [
        ([(2, (1, 2, 4), (0, 0, 0)), (3, (1, 1), (0, 1))], 1, 0),
        ([(5, (1, 2), (3, 1)), (3, (1, 2), (0, 1))], 7, 3),
        ([(2, (3, 2), (0, 1)), (5, (3, 2), (0, 1))], 3, 1),
    ]
    found = 0
    for params, modulus, residue in cases:
        graphs = [CarryGraph(*g) for g in params]
        expected = []
        for t in range(2, 1501):
            if t % modulus != residue:
                continue
            for p, slopes, offsets in params:
                entries = [a*t+b for a, b in zip(slopes, offsets)]
                left = sum(entries)
                coefficient = 1
                for x in entries:
                    coefficient *= comb(left, x)
                    left -= x
                if coefficient % p == 0:
                    break
            else:
                expected.append(t)
        actual, _, _ = intersection(graphs, 1500, modulus, residue)
        assert actual == expected
        found += len(actual)
    assert found > 0
    print(f"Direct multinomial comparison agrees on {found} accepted inputs across three tests.")


def dimension_percent(g):
    length = 64
    bound = g.block_bound(length)**100
    return next(j for j in range(120) if bound <= g.base**(j*length))


def factor_five_seven(limit):
    for d, part, modulus, residue in [(5, 3, 15, 7), (7, 6, 72, 11)]:
        graphs = [CarryGraph(p, (part, d-part), (0, 1)) for p in (2, 3, d)]
        actual, nodes, tested = intersection(graphs, limit, modulus, residue)
        assert not actual, (d, actual)
        print(f"Factor {d}: {nodes} intervals, {tested} direct checks, zero survivors.", flush=True)
        if d == 5:
            h = 16384
            b2 = graphs[0].block_bound(h).bit_length()
            b5 = graphs[2].block_bound(h).bit_length()
            assert (b2, b5) == (9789, 19139) and 5**1000 > 2**2321
            assert 1000*(2321*b2+1000*b5) < 1101*2321*h

    # At odd weights with full factor 7, the row covering the predecessor
    # also avoids 2 and the prime scale u. Its part is k*(N-1)/7.
    for u, k in ((3, 3), (3, 6), (5, 5)):
        modulus = 14*u
        residue, = [t for t in range(0, modulus, 2)
                    if (7*t+1) % u == 0 and k*t % 7 == 1]
        graphs = [CarryGraph(p, (k, 7-k), (0, 1)) for p in (2, u, 7)]
        actual, nodes, tested = intersection(graphs, limit, modulus, residue)
        expected = [8452] if k == 5 and limit >= 8452 else []
        assert actual == expected, (u, k, actual)
        for t in actual:
            m = 7*t
            parts = [(m//w)*pow(m//w, -1, w) for _, w in prime_powers(m)]
            assert len(parts) == 3 and sum(parts) != m+1
        if k in (3, 5):
            h = 16384
            b2 = graphs[0].block_bound(h).bit_length()
            b7 = graphs[2].block_bound(h).bit_length()
            assert 7**1000 > 2**2807
            assert 1000*(2807*b2+1000*b7) < 1212*2807*h
        print(f"Odd factor 7, scale {u}, multiplier {k}: {nodes} intervals, "
              f"{tested} direct checks, aggregate survivors {actual}.", flush=True)


def bounded_cofactors(limit):
    spf = smallest_factors(1000)
    total = nodes_all = tested_all = 0
    worst = 0
    for d in COFACTORS:
        factors = atoms(d, spf)
        for shape in shapes(d, factors):
            if not survives(d, factors, shape, spf, False):
                continue
            u, a, b, c, *_ = shape
            primes = forced_bases(shape, factors, atoms(u, spf))
            graphs = [CarryGraph(p, (a, b, c), (0, 0, 1)) for p in primes]
            bounds = sorted((dimension_percent(g), g.base) for g in graphs)
            (e, p), (f, q) = bounds[:2]
            assert p != q and e+f < 100
            worst = max(worst, e+f)
            actual, nodes, tested = intersection(graphs, limit, u, -pow(d, -1, u))
            assert not actual, (d, shape, actual)
            total += 1
            nodes_all += nodes
            tested_all += tested
            print(f"d={d} slopes={(a,b,c)}: bases {p},{q} at most {e}+{f} percent, "
                  f"{nodes} intervals, {tested} direct checks.", flush=True)
    assert total == 16
    print(f"{total} shapes, largest dimension sum {worst}/100, {nodes_all} intervals, "
          f"{tested_all} direct checks, zero survivors.")


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--digits", type=int, default=1000,
                        help="inclusive upper limit is ten to this power")
    args = parser.parse_args()
    if args.digits < 1:
        parser.error("digits must be positive")
    small_crosscheck()
    factor_five_seven(10**args.digits)
    bounded_cofactors(10**args.digits)
