"""Check the one four-prime list for all 2**76 specified weights."""
from itertools import combinations_with_replacement
from math import comb, gcd

from check_cluster_certificates import is_prime
from check_cluster_minima import candidate_bases, composite_witness
from check_quartic_triangle import digits


def splits(total):
    for a in range(total+1):
        for b in range(total-a+1):
            yield a, b, total-a-b


def low_block_controls():
    rows = forced = joint = 0
    for w, p in ((4, 2), (5, 5), (7, 7), (8, 2), (9, 3), (11, 11)):
        for f in range(2*w, 101, w):
            for gap in range(1, 5):
                r = f-gap
                if not is_prime(r) or r == p:
                    continue
                for d in range(1, 5):
                    if d*gap+1 >= w or d >= r:
                        continue
                    for epsilon in (0, 1):
                        n = d*f+epsilon
                        for upper in splits(d):
                            for lower in splits(d*gap+epsilon):
                                parts = tuple(e*r+j for e, j in zip(upper, lower))
                                if not min(parts):
                                    continue
                                a, b, _ = parts
                                coef = comb(n, a)*comb(n-a, b)
                                assert coef % r
                                rows += 1
                                residues = tuple(x % w for x in parts)
                                if sum(residues) != epsilon:
                                    assert coef % p == 0
                                    continue
                                forced += 1
                                assert parts == tuple(e*f+h for e, h in zip(upper, residues))
                                assert sum(e > 0 for e in upper) >= 2
                                joint += gcd(coef, p*r) == 1
    assert 0 < joint <= forced < rows
    print(f"{rows} exact small coefficients, {forced} forced row forms, {joint} joint avoidances.")


def example():
    p = 2**256
    a, b = 297, 3072559256396398118163687
    q, r = p+a, p+b
    f = p*q*r
    last = 2**75
    v = f-1185
    assert all(is_prime(x) for x in (q, r, v))
    assert 0 < 2*a < b and q*r % last == last-1
    assert 2*last*b*b <= p and 1185*last+1 < p and last < v
    assert 180*last*b*b <= p  # Independent three-power criterion.
    assert len({2, q, r, v}) == 4
    for e, h in combinations_with_replacement((1, 2, 3, last//2, last-2, last-1), 2):
        if e+h <= last:
            assert e*q*r % last + h*q*r % last >= last
            assert e*p*r % q + h*p*r % q > q
            assert (e*p*q//r) % r + (h*p*q//r) % r > r
    value = 1
    for _ in range((last*f).bit_length()):
        value = 2*value % q
        assert value != q-1
    assert (q+1) % p and (r+1) % p and last < p-1
    for n in (f, f+1):
        assert all(composite_witness(base) for base in set(candidate_bases(n)))
        assert sum(digits(n, 2)) > 2
    assert gcd(f, 225) == 1
    def exactly(n, prime):
        return n % prime == 0 and n % (prime*prime) != 0
    residues = [d for d in range(225) if
                (exactly(d*f-1, 3) and exactly(d*f+1, 5)) or
                (exactly(d*f-1, 5) and exactly(d*f+1, 3))]
    assert len(residues) == 16 and 0 not in residues
    minimal_pairs = sum(1+(last-d)//225 for d in residues)
    print(f"One list certifies {2*last} weights, with {minimal_pairs} pairs minimal by the wheel.")
    print(f"At d=1 both minimum sizes are also four. The fourth prime is F-1185.")


if __name__ == "__main__":
    low_block_controls()
    example()
