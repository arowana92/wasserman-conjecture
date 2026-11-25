"""Verify the factor-five binary enumeration independently and check exact valuations."""

from argparse import ArgumentParser
from fractions import Fraction
from pathlib import Path
from subprocess import run


def generated(executable, bits, max_k):
    output = run([executable, str(bits), str(max_k)], capture_output=True,
                 text=True, check=True)
    triples = [tuple(map(int, line.split())) for line in output.stdout.splitlines()]
    assert len(set(triples)) == len(triples)
    return set(triples), output.stderr.strip()


def independent_small_box():
    triples = set()
    for q in range(7, 2**16, 2):
        if q % 3 == 0 or q % 5 == 0 or ((3*q) & (2*q)):
            continue
        for k in range(3, 513):
            p = 1 << k
            n = 5*p*q+1
            if n % 3 or n % 9 == 0:
                continue
            l = pow(5*p, -1, q)
            row = (3*p*q, 5*p*l, n-3*p*q-5*p*l)
            if min(row) < 1 or any(x % 3 for x in row):
                continue
            if [x % 5 for x in row] != [1, 0, 0]:
                continue
            a, b, c = row
            if (a & b) | (a & c) | (b & c):
                continue
            triples.add((q, k, l))
    return triples


def valuation(row, p):
    n = sum(row)
    power = p
    exponent = first = 0
    while power <= n:
        contribution = n//power - sum(x//power for x in row)
        if contribution and not first:
            first = power
        exponent += contribution
        power *= p
    return exponent, first


def main():
    parser = ArgumentParser()
    parser.add_argument("generator", type=Path)
    executable = str(parser.parse_args().generator.resolve())
    small = independent_small_box()
    assert small and small == generated(executable, 16, 512)[0]

    triples, summary = generated(executable, 42, 2048)
    assert len(triples) == 107
    maximum = {3: 0, 5: 0}
    for q, k, l in triples:
        p = 1 << k
        n = 5*p*q+1
        a, b = 3*p*q, 5*p*l
        row = (a, b, n-a-b)
        assert min(row) > 0 and n % 3 == 0 and n % 9 != 0
        assert not any(x % 3 for x in row)
        assert [x % q for x in row] == [0, 1, 0]
        assert [x % 5 for x in row] == [1, 0, 0]
        assert not ((a & b) | (a & row[2]) | (b & row[2]))
        for base in maximum:
            value, first = valuation(row, base)
            assert value > 0
            maximum[base] = max(first, maximum[base])
    assert maximum == {3: 3**8, 5: 5**4}
    rho = Fraction(303, 200)
    assert rho**4-rho**3-rho**2+rho-1 > 0 and rho**5 < 8
    print(f"Independent enumeration agrees on {len(small)} small-box candidates.")
    print(summary)
    print(f"All {len(triples)} rows carry in both bases. Maximum first layers: {maximum}.")


if __name__ == "__main__":
    main()
