"""Rule out perfect-power cofactors by successively lifting six exponents."""
from itertools import combinations_with_replacement, groupby, product
from math import gcd

from flint import fmpz, fmpz_mod_ctx


# A prime may be used from the first period h at which 7^h = 1 modulo it.
# The lists need not contain every possible prime divisor of 7^h-1.
PRIMES = {
    5: (11, 31, 61, 191, 281, 541, 2801, 4021, 9901, 154351, 159871, 295831, 6568801, 7777981,
         555915824341),
    7: (29, 113, 911, 4733, 51031, 309079, 11898664849),
    11: (23, 67, 661, 1123, 1409, 3631, 29569, 293459, 1532917, 2784937, 10746341, 12323587,
         83960385389, 1357105535093947),
    13: (53, 79, 486643, 7524739, 23850061, 44975113, 228511817, 16148168401, 115868130379),
    17: (103, 14009, 365773, 1135261, 12192877, 2767631689, 29078814248401, 90920746649833),
    19: (419, 457, 19609, 351121, 879399649, 4058036683, 24006520939, 4534166740403,
         6957533874046531, 275430251632699963),
    23: (47, 139, 277, 3083, 26681, 649981, 15513961, 31479823396757, 3421093417510114543),
    29: (59, 127540261, 2576743207, 5406066031, 71316922984999, 196915704073465747,
         1284564941030772067, 34720058992462423813, 358475907408445923469,
         402488219476647465854701),
    31: (311, 373, 21143, 149731, 5463751, 9754399, 5420506947192709, 3999088279399464409),
    37: (149, 223, 1999, 2887, 6217, 33301, 177379, 28063132333, 4805345109492315767981401,
         15572244900182528777225808449),
    41: (83, 267403, 4849399, 8999993, 20515909, 395617447, 1636258751, 78009515593,
         4362139336229068656094783),
    43: (173, 947, 1033, 10837, 14449, 19867, 20899, 163573, 235211, 3799051, 4084141, 4749781,
         384439093, 166003607842448777, 2192537062271178641),
    47: (283, 7333, 5960917, 13722816749522711, 63681511996418550459487,
         655417289594537954307682339857743931943,
         285917932930729262773854471418794183208735026049402141999388275408759),
}
STEPS = (1, 2, 6, 12, 36)


def lifts(values, old_h, h):
    groups = []
    for r, xs in groupby(values):
        multiplicity = sum(1 for _ in xs)
        groups.append(combinations_with_replacement(range(r, h, old_h), multiplicity))
    for blocks in product(*groups):
        yield tuple(x for block in blocks for x in block)


def check(ell):
    available = sorted(PRIMES[ell])
    assert all(fmpz(r).is_prime() and (r-1) % ell == 0 for r in available)
    live = None
    old_h = 1
    counts = []
    for step in STEPS:
        h = step*ell
        active = [r for r in available if pow(7, h, r) == 1]
        available = [r for r in available if r not in active]
        tables = {r: {pow(x, ell, r) for x in range(r)} for r in active if r < 200000}
        rings = {r: fmpz_mod_ctx(r) for r in active if r not in tables}
        # Each table tests the cofactor t directly, including zero.
        mod = gcd(432, 7**h-1)//6
        repunits = [(7**j-1)//6 for j in range(h)]
        candidates = (combinations_with_replacement(range(h), 6) if live is None else
                      (xs for old in live for xs in lifts(old, old_h, h)))
        new = []
        for xs in candidates:
            t = sum(repunits[j] for j in xs)
            if t % mod != 11 % mod:
                continue
            for r in active:
                x = t % r
                if r in tables:
                    if x not in tables[r]:
                        break
                elif x and rings[r](x)**((r-1)//ell) != 1:
                    break
            else:
                new.append(tuple(sorted(xs)))
        counts.append(len(new))
        if not new:
            print(ell, counts, flush=True)
            return
        live, old_h = new, h
    raise ValueError(("periods exhausted", ell, counts))


if __name__ == "__main__":
    from sys import argv
    for exponent in map(int, argv[1:]) if argv[1:] else PRIMES:
        check(exponent)
