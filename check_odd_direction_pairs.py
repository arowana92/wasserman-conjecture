"""Run the parity enumeration and verify its records by two exact valuation formulas."""

from argparse import ArgumentParser
from concurrent.futures import ThreadPoolExecutor
from hashlib import sha256
from math import gcd, isqrt, prod
from pathlib import Path
from subprocess import check_output

from check_direction_carry_count import no_carry
from check_direction_pair_count import enumerate_pairs, valuation
from check_primitive_row_scales import direct_layers, powers, profile, scaled


def records_from(text):
    result = set()
    for line in text.splitlines():
        n, a, b, c, d, e, f = map(int, line.split())
        result.add((n, (a, b, c), (d, e, f)))
    return result


def native_records(binary, limit, jobs):
    if jobs == 1:
        return records_from(check_output([binary, str(limit)], text=True))

    def shard(d):
        return records_from(check_output([binary, str(limit), str(d), str(d+3)], text=True))

    with ThreadPoolExecutor(max_workers=jobs) as pool:
        return set().union(*pool.map(shard, range(2, isqrt(limit//3)+1, 4)))


def canonical_text(records):
    return "".join(" ".join(map(str, (n, *x, *y)))+"\n" for n, x, y in sorted(records))


def digit_sum(z, p):
    total = 0
    while z:
        z, digit = divmod(z, p)
        total += digit
    return total


def minimum_cover(edges):
    # Any cover must hit a smallest remaining edge. A disjoint subfamily
    # lower-bounds the remaining size in each branch.
    reduced = []
    for edge in sorted(set(edges), key=lambda e: (len(e), sorted(e))):
        if not any(s <= edge for s in reduced):
            reduced.append(edge)

    def visit(rest, slots):
        if not rest:
            return ()
        used = set()
        lower = 0
        for edge in rest:
            if used.isdisjoint(edge):
                lower += 1
                used.update(edge)
        if lower > slots:
            return None
        for p in sorted(rest[0], key=lambda p: (-sum(p in e for e in rest), p)):
            tail = visit([e for e in rest if p not in e], slots-1)
            if tail is not None:
                return (p, *tail)
        return None

    for k in range(len(reduced)+1):
        cover = visit(reduced, k)
        if cover is not None:
            return tuple(sorted(cover))
    raise ValueError("An empty prime set has no transversal")


def exclude_covers_without_two(limit):
    # The six endpoint coefficients below have disjoint odd prime supports.
    witnesses = [
        ((212, 464, 1537), (491, 631, 631), 3879389),
        ((316, 456, 1501), (551, 859, 859), 5157437),
        ((3, 743, 743), (68, 1168, 1241), 3688253),
        ((95, 131, 393), (104, 824, 1339), 1403273),
        ((44, 604, 1661), (315, 359, 359), 2385197),
        ((184, 344, 989), (375, 617, 617), 2440853),
    ]
    seen = set()
    count = 0
    for x, y, n in witnesses:
        if max(sum(x), sum(y)) > limit:
            continue
        assert profile(x, y)[2] == n and gcd(*x) == gcd(*y) == 1
        odd = {p for p, _ in powers(n)+powers(n-1) if p > 2}
        assert odd and seen.isdisjoint(odd)
        seen.update(odd)
        count += 1
    return count


def report(records, limit):
    cases = []
    for n, x, y in records:
        assert x <= y and gcd(*x) == gcd(*y) == 1 and max(sum(x), sum(y)) <= limit
        assert profile(x, y)[2] == n
        rows = scaled(x, n), scaled(y, n)
        assert direct_layers(n, rows) and any(no_carry(2, row) for row in rows)
        factors = {}
        for p, w in powers(n)+powers(n-1):
            exponents = [valuation(n, row, p) for row in rows]
            for e, row in zip(exponents, rows):
                assert e*(p-1) == sum(digit_sum(a, p) for a in row)-digit_sum(n, p)
            if min(exponents):
                factors[p] = min(w, p**min(exponents))
        if not factors:
            raise RuntimeError(f"Exact valuation candidate with endpoint gcd one: {n}, {x}, {y}")
        cases.append((max(sum(x), sum(y)), n, x, y, prod(factors.values()), frozenset(factors)))
    cases.sort()
    cover = (2, *minimum_cover([r[-1] for r in cases]))
    disjoint = exclude_covers_without_two(limit)
    if limit == 160000:
        assert len(cases) == 272 and len({r[1] for r in cases}) == 195
        assert cover == (2, 3, 5, 7, 11, 23) and disjoint == 6
    print(f"Primitive sums through {limit}: {len(cases)} parity survivors at "
          f"{len({r[1] for r in cases})} distinct weights.")
    print("Smallest cover containing 2:", cover)
    print("Lower bound if 2 is omitted:", disjoint)
    if cases:
        print("Smallest gcd among parity survivors:", min(r[4] for r in cases))
        print("Largest parity-survivor weight:", max(r[1] for r in cases))
        for primes in ((2,), (2,3), (2,5), (2,3,5), (2,3,5,7), (2,3,5,7,11)):
            first = next((r for r in cases if r[-1].isdisjoint(primes)), None)
            print("First missed by", primes, ":", None if first is None else first[:-1])
        heights = sorted({r[0] for r in cases})
        for k in range(2, len(cover)+1):
            lo, hi = 0, len(heights)-1
            while lo < hi:
                mid = (lo+hi)//2
                subcover = minimum_cover([r[-1] for r in cases if r[0] <= heights[mid]])
                if 1+len(subcover) >= k:
                    hi = mid
                else:
                    lo = mid+1
            print(f"First sum bound requiring {k} primes in a cover containing 2: {heights[lo]}")
    print("Canonical SHA-256:", sha256(canonical_text(records).encode("ascii")).hexdigest())
    return cover


def compare(records, limit):
    cropped = {r for r in records if max(sum(r[1]), sum(r[2])) <= limit}
    all_pairs = enumerate_pairs(limit)[1]
    if limit >= 149:
        assert min(max(sum(x), sum(y)) for x, y, n in all_pairs) == 149
    expected = {(n, x, y) for x, y, n in all_pairs
                if any(valuation(n, scaled(z, n), 2) == 0 for z in (x,y))}
    assert expected == cropped
    print(f"Complete Python enumeration agrees on all {len(cropped)} records through {limit}.")


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("--limit", type=int, default=6000)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--compare", type=int, help="also enumerate every lowest-layer pair up to this sum")
    parser.add_argument("--output", type=Path, help="save sorted records")
    args = parser.parse_args()
    if not 3 <= args.limit <= 160000 or args.jobs < 1:
        parser.error("the sum bound must be 3..160000 and jobs must be positive")
    if args.compare is not None and not 3 <= args.compare <= args.limit:
        parser.error("the comparison bound must be between 3 and the sum bound")
    records = native_records(str(args.binary.resolve()), args.limit, args.jobs)
    report(records, args.limit)
    if args.compare is not None:
        compare(records, args.compare)
    if args.output:
        args.output.write_text(canonical_text(records), encoding="ascii")
