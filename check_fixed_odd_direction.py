"""Check the fixed binary direction construction, including mixed and repeated prime factors."""

from math import gcd, prod

from check_direction_carry_count import no_carry
from check_direction_pair_count import valuation
from check_prescribed_common_prime import crt, fill, is_prime, word_for
from check_primitive_row_scales import direct_layers, powers, profile


def data(u, bound):
    primes = [p for p in range(3, bound+1, 2) if is_prime(p)]
    excluded = [p for p in primes if u % p]
    base, period = crt([(1, 2), *((u//w * pow(2, -1, w), w) for _, w in powers(u))])
    omit = prod(excluded)
    d = base+period*max(0, (u-base)//period+1)
    while gcd(d, omit) != 1:
        d += period
    e = d+period
    while gcd(d*omit, e) != 1:
        e += period
    return d, e, excluded


def target(u, d, e, excluded):
    D = d*e
    r0, _ = crt([(pow(u*e, -1, d), d), (pow(u*d, -1, e), e)])
    return crt([(1, u), ((1+2*D*r0)*pow(u, -1, D*D), D*D),
                *((-pow(u, -1, p), p) for p in excluded)])


def verify(x, d, e, bound, v):
    u = sum(x)
    n = u*v
    r, rem = divmod(n-1, 2*d*e)
    assert not rem and r>0 and gcd(r, d*e) == gcd(d, e) == 1
    first = tuple(v*a for a in x)
    second = (u*e*r, u*d*r, (2*d*e-u*(d+e))*r+1)
    y = tuple(a//u for a in second)
    assert gcd(*first) == v and gcd(*second) == u and sum(y) == v and min(y) > 0
    assert gcd(u, v) == 1 and profile(x, y)[2] == n and no_carry(2, first)
    assert [sum(a % m for a in second) for m in (d, e, r)] == [1, 1, 1]
    part = 1
    for p in range(2, bound+1):
        if not is_prime(p):
            continue
        power, endpoint = 1, n*(n-1)
        while endpoint % p == 0:
            power *= p
            endpoint //= p
        exponent = min(valuation(n, first, p), valuation(n, second, p))
        part *= min(power, p**exponent)
    assert part == u
    return n, first, second


def prescribed_examples():
    for x, bound, v in (
        ((1, 4, 16), 7, 1462974881671711613292555143881767096792419140609),
        ((1, 2, 24), 3, 98080088852629449863035390360723719343166218872253841473),
    ):
        u = sum(x)
        d, e, excluded = data(u, bound)
        r0, modulus = target(u, d, e, excluded)
        assert v % modulus == r0
        n, _, second = verify(x, d, e, bound, v)
        for p, w in powers(u):
            assert [(a//w) % w for a in second[:2]] == [w-1, w-1]
            assert 2*d % w == 2*e % w == (u//w) % w
        print(f"Primitive sum {u}: d={d}, e={e}, v={v}, weight={n}.")

    # Smaller scalar for the (1,2,4) example. Its base-seven carry is direct.
    u, d, e = 7, 11, 13
    residue, modulus = target(u, d, e, [3, 5])
    layers = fill(modulus, 3, 50)
    v = 1+8*word_for(layers, modulus, 3, (residue-1)*pow(8, -1, modulus) % modulus)
    n, first, second = verify((1, 2, 4), d, e, 7, v)
    assert n == 15874320294119
    common = prod(min(w, p**min(valuation(n, row, p) for row in (first, second)))
                  for p, w in powers(n)+powers(n-1))
    assert common == 16319685886490987
    print(f"Primitive sum 7: v={v}, weight={n}, exact endpoint gcd={common}.")


def repair_controls():
    tested = 0
    for x in ((1, 2, 4), (1, 4, 16), (5, 8, 18)):
        gap = sum(x).bit_length()
        high = gap+12
        for modulus in (3, 7, 15, 31):
            available = [1 << (high+j*gap) for j in range(modulus-1)]
            words = {0: 0}
            for bit in available:
                old = list(words.items())
                for residue, word in old:
                    words.setdefault((residue+bit) % modulus, word+bit)
                assert len(words) > len(old) or len(words) == modulus
            assert len(words) == modulus
            for w in range(128):
                if not no_carry(2, tuple(a*w for a in x)):
                    continue
                bottom = 1+(w << gap)
                for r in range(modulus):
                    v = bottom+words[(r-bottom) % modulus]
                    assert v % modulus == r and (v % (1 << high)) >> gap == w
                    assert no_carry(2, tuple(a*v for a in x))
                    tested += 1
    print(f"All {tested} block repairs preserved the free word and disjoint digits.")


def arbitrary_directions():
    tests = 0
    for x in ((1, 1, 1), (1, 1, 2), (1, 2, 3), (1, 2, 4), (2, 3, 4), (3, 5, 8)):
        u = sum(x)
        d = 2*u+1
        e = 2*u*d+1
        D = d*e
        r0, _ = crt([(pow(u*e, -1, d), d), (pow(u*d, -1, e), e)])
        start, step = crt([(1, u), ((1+D*r0)*pow(u, -1, D*D), D*D)])
        for j in range(4):
            v = start+j*step
            n = u*v
            r = (n-1)//D
            row = (u*e*r, u*d*r, (D-u*(d+e))*r+1)
            y = tuple(a//u for a in row)
            assert gcd(*row) == u and sum(y) == v and min(y) > 0
            assert profile(x, y)[2] == n and direct_layers(n, (tuple(v*a for a in x), row))
            tests += 1
    print(f"Checked {tests} assignments for arbitrary fixed directions.")


if __name__ == "__main__":
    prescribed_examples()
    repair_controls()
    arbitrary_directions()
