"""Independent interval cover by Legendre factorial floors and Arb logarithms."""
from argparse import ArgumentParser
from functools import lru_cache

from flint import arb, ctx
from check_logarithmic_intervals import forced_configurations

BITS = 640
UNIT = 1 << BITS


@lru_cache(maxsize=None)
def logarithm(p):
    return arb(p).log()


def layers(p, j, max_slope, margin_bits=None):
    """Outward intervals for consecutive powers of p in units UNIT / 2**j."""
    if p == 2:
        lo = hi = 2*UNIT
        exponent = j+1
    else:
        lp = logarithm(p)
        l2 = logarithm(2)
        k = (((j+1)*l2/lp).floor()).unique_fmpz()
        assert k is not None
        exponent = int(k) + 1
        value = (exponent*lp - j*l2).exp() * UNIT
        lo = int(value.lower().floor().unique_fmpz())
        hi = int(value.upper().ceil().unique_fmpz())
        assert hi-lo <= 2
    while lo < 2*UNIT*max_slope:
        lo, hi, exponent = p*lo, p*hi, exponent+1
    if margin_bits is not None and p != 2:
        radius = -(-hi // (1 << margin_bits))
        lo, hi = lo-radius, hi+radius
    result = []
    while lo > 1 and exponent >= 1:
        result.append((lo, hi))
        lo, hi, exponent = lo//p, -(-hi//p), exponent-1
    return result


def factorial_cover(slopes, offsets, tables, j):
    """Cover closed real intervals by strict positive factorial-floor differences."""
    # In normalized coordinates an offset of one is UNIT/2**j.
    unit_offset = 1 << max(0, BITS-j)
    coeffs = tuple(zip(slopes, (b*unit_offset for b in offsets)))
    largest = max(slopes)
    total = sum(slopes)
    left = UNIT
    steps = 0
    while left < 2*UNIT:
        best = left
        for table in tables:
            for lo, hi in table:
                # A single floor for the largest part persists over width at most hi/largest.
                if hi <= largest*(best-left):
                    break
                minimum = (total*left)//hi
                if not minimum:
                    continue
                floors = [(a*left+b)//lo for a,b in coeffs]
                if minimum > sum(floors):
                    # For x <= right the upper part estimates stay strictly
                    # below the following multiples of this layer's lower bound.
                    right = min(((f+1)*lo-b-1)//a for f,(a,b) in zip(floors,coeffs))
                    best = max(best, right)
        if best <= left:
            return False, steps
        left = best
        steps += 1
    return True, steps


def controls():
    from math import comb
    for j in (1, 7, 80, 300, 512):
        for p in (2, 3, 7, 19):
            prev = None
            for lo, hi in layers(p, j, 10):
                if prev is None:
                    power = 1
                    while power*UNIT < lo*(1<<j):
                        power *= p
                else:
                    power = prev//p
                assert lo*(1<<j) <= power*UNIT <= hi*(1<<j)
                prev = power

    examples = empty = 0
    for j in range(7):
        for slopes, offsets, ps in (((3,2),(0,1),(2,5)), ((1,1),(0,1),(2,3))):
            ok, _ = factorial_cover(slopes, offsets, [layers(p,j,sum(slopes)) for p in ps], j)
            accepted = sum(all(comb(sum(slopes)*t+1, slopes[0]*t) % p for p in ps)
                           for t in range(1<<j, 1<<(j+1)))
            assert not (ok and accepted)
            examples += accepted
            empty += ok
    assert examples and empty
    print(f"Direct controls: {examples} positive cases and {empty} certified empty bins.", flush=True)


def run(start, count, include_cuts):
    ctx.prec = start.bit_length() + count.bit_length() + BITS + 120
    controls()
    rows = {(sl,(0,0,1),ps) for _,sl,ps in forced_configurations()}
    if include_cuts:
        for d,k,u in ((5,3,3),(7,3,3),(7,5,5),(7,6,3)):
            rows.add(((k,d-k),(0,1),(2,u,d)))
    widths = {p: max(sum(sl) for sl,_,ps in rows if p in ps)
              for _,_,ps in rows for p in ps}
    overall = worst = 0
    for j in range(start, start+count):
        tables = {p: layers(p,j,width,440) for p,width in widths.items()}
        for sl,off,ps in rows:
            ok, steps = factorial_cover(sl, off, [tables[p] for p in ps], j)
            assert ok, (j,sl,ps)
            overall += steps
            worst = max(worst, steps)
    print(f"Certified {count} bins, {len(rows)} rows each, {overall} factorial-floor intervals,"
          f" maximum {worst}.", flush=True)


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--start', type=int, default=10**100)
    parser.add_argument('--count', type=int, default=16384)
    parser.add_argument('--include-cuts', action='store_true')
    args = parser.parse_args()
    run(args.start,args.count,args.include_cuts)
