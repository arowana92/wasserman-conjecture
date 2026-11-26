"""Independent finite check using exponent directions and direct integer gcds."""
from math import gcd, isqrt


def residual_masks(g, budget, denominator, slope_cost):
    """Uncovered f modulo g for each e. The cost uses integer numerator units."""
    all_bits = (1 << g)-1
    uncovered = [all_bits]*g
    for n in range(1, min(g-1, budget//slope_cost)+1):
        width = (budget-n*slope_cost)//denominator
        bins = [0]*g
        for f in range(g):
            bins[n*f % g] |= 1 << f
        prefix = [0]
        for row in bins:
            prefix.append(prefix[-1] | row)
        for e in range(g):
            re = n*e % g
            if re > g//2:
                re -= g
            if abs(re) > width:
                continue
            low, high = max(0,re)-width, min(0,re)+width
            # Permitted signed residues range from low<=0 to high>=0.
            covered = prefix[high+1]
            if low < 0:
                covered |= all_bits ^ prefix[g+low]
            uncovered[e] &= ~covered
        if not any(uncovered):
            break
    return uncovered


def run(binary):
    bases = (2,) if binary else (3,5,7,11,13,17,19,23)
    total = direct = large = 0
    for g in range(68 if binary else 26, 1383 if binary else 900):
        if binary:
            if g & (g-1) and any(g % k == 0 for k in range(2,isqrt(g)+1)):
                continue
            residual = residual_masks(g, 5*g, 30, 48)
            m = (3*g-1)//2
        else:
            if g % 3 == 0:
                continue
            exponent = 4*max(101,3*g)-9*g
            residual = residual_masks(g, exponent-18, 18, 18)
            m = 2*g-1
        if not any(residual):
            continue
        for p in bases:
            powers = [1]
            for _ in range(m):
                powers.append(p*powers[-1])
            modulus = 1+powers[g]
            modulus //= modulus & -modulus
            for e in range(int(binary), m+1):
                start = e+int(binary)
                stop = m
                if not binary:
                    # Rotate the three residues modulo 2g past a largest gap.
                    # That last gap 2g-f is at least both e and f-e.
                    stop = min(stop, 2*g-e, (2*g+e)//2)
                if start > stop:
                    continue
                total += stop-start+1
                mask = residual[e%g]
                mask = ((mask | (mask << g)) >> start) & ((1 << (stop-start+1))-1)
                f = start-1
                while mask:
                    jump = (mask & -mask).bit_length()
                    f += jump
                    mask >>= jump
                    common = gcd(modulus, 1+powers[e]+powers[f])
                    direct += 1
                    if binary:
                        h = max(101,g,f)
                        assert common <= powers[(4*h-3*g)//6]
                    elif common > powers[max(0, exponent//18-1)]:
                        large += 1
                        assert 256**6*common**18 < 27**6*p**exponent
    print(f"{bases}: {direct} direct numerical gcds, {large} above the rounded cutoff", flush=True)
    assert (direct,large) == ((4091854,0) if binary else (22297376,1))
    assert total > direct


if __name__ == "__main__":
    run(True)
    run(False)
