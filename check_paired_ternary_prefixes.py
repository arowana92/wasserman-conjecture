"""Independent reverse enumeration at remainder 11. Requires python-flint."""

from bisect import bisect_left
import re

from flint import fmpz

from check_paired_cofactor_range import paired_weights
from check_residual_power_classes import paired_binary_ternary_cofactors


def ternary_suffixes(digits, powers):
    # A qualifying E beginning at ternary position h has no nonzero digits
    # in positions h+1 through h+4 of the inverse, since 121 > 3^4.
    # Reconstruct larger suffixes by appending only the newly exposed digits.
    digits = '0000' + digits
    suffix = fmpz(0)
    used = 0
    for match in reversed(list(re.finditer(r'0{4}(?=[12])', digits))):
        start = match.end()
        prefix = digits[start:len(digits)-used]
        # Convert in short pieces also when a run exceeds Python's default
        # single-string conversion limit.
        extension = 0
        for index in range(0, len(prefix), 1000):
            piece = prefix[index:index+1000]
            extension = extension*3**len(piece) + int(piece,3)
        suffix += extension*powers[used]
        used = len(digits)-start
        yield suffix


def reverse_enumeration(limit, classes):
    limit = fmpz(limit)
    powers = [fmpz(1)]
    while powers[-1] < 121*limit:
        powers.append(3*powers[-1])
    modulus = powers[-1]
    max_b = bisect_left(powers, limit)

    # 3 has exact order 2^(k-2) modulo 2^k for k >= 3.
    last_k = 3
    while 2**(last_k-2) < max_b:
        last_k += 1
    logs = {}
    for k in range(3,last_k+1):
        m = 2**k
        residue = 1
        for b in range(1,min(max_b,2**(k-2))+1):
            residue = 3*residue % m
            logs[k,residue] = b

    inverse = fmpz(1)
    P = fmpz(1)
    rows = set()
    for a in range(1,limit.bit_length()+1):
        P *= 2
        upper = (limit-1)//(33*P)
        if upper < 1:
            break
        inverse = (inverse + (modulus if inverse % 2 else 0))//2
        end = bisect_left(powers,121*upper)
        digits = (inverse % powers[end]).str(base=3)
        for E in ternary_suffixes(digits,powers):
            if E > upper:
                break
            if int(E % 1320) not in classes:
                continue
            level = powers[bisect_left(powers,121*E)]
            if inverse % level != E:
                continue
            binary_exp = (121*E-1).bit_length()
            k = min(last_k,binary_exp)
            residue = pow(int(E % (2**k)),-1,2**k)
            first_b = logs.get((k,residue))
            if first_b is None:
                continue
            for b in range(first_b,max_b+1,2**(k-2)):
                if 1+11*E*P*powers[b] > limit:
                    break
                if (E*powers[b]-1) % (fmpz(2)**binary_exp) == 0:
                    rows.add((int(E),a,b))
    return rows


if __name__ == '__main__':
    limit = 10**10000
    classes = paired_binary_ternary_cofactors()
    rows = reverse_enumeration(limit,classes)
    assert rows == paired_weights(2,3,11,limit,classes)
    assert len(rows) == 71
    print('The binary and ternary inverse searches agree on all 71 tuples.')
