"""Cross-check fixed binary exponents with a factorial-layer automaton."""

import argparse
from collections import deque
from pathlib import Path
import subprocess


def templates(k):
    modulus = 2 << k
    prefixes = [1]
    bit = 2
    while bit < modulus:
        prefixes = [n for v in prefixes for n in (v, v + bit)
                    if not (3*n & 2*n & bit)]
        bit *= 2
    result = {}
    for v in prefixes:
        c = pow(v, -1, modulus)
        if c % 15 == 0:
            result[c] = v
    return result


def carry_free(total, first, second, p):
    # The first two parts fit into the total digit by digit exactly when
    # the third (their complement) adds to them without any carry.
    if p == 2:
        return not ((first & second) or ((first | second) & ~total))
    while total or first or second:
        if first % p + second % p > total % p:
            return False
        total, first, second = total//p, first//p, second//p
    return True


def graph(k, c, q, p, whole=False):
    P = 1 << k
    if p == 2:
        m, residues = 2*P, [q]
    elif p == 3:
        m = 9
        residues = [v for v in range(m) if (5*P*v+1) % 3 == 0 and (5*P*v+1) % 9]
    else:
        m, residues = 5, [2 * pow(P, -1, 5) % 5]
    pending = []
    for v in residues:
        n, a, b = 5*P*v+1, 3*P*v, c*v
        if carry_free(n % m, a % m, b % m, p):
            state = (n//m, a//m, b//m, (3*v)//m if p==2 else 0, (2*v)//m if p==2 else 0)
            pending.append(state)
    start = pending.copy()
    seen = set(pending)
    reverse = {}
    terminal = []
    for state in pending:
        n, a, b, f, g = state
        if carry_free(n, a, b, p) and not (f & g):
            if not whole:
                return True
            terminal.append(state)
        for digit in range(p):
            nn, aa, bb = n+5*P*digit, a+3*P*digit, b+c*digit
            if aa % p + bb % p > nn % p:
                continue
            ff, gg = (f+3*digit, g+2*digit) if p==2 else (0, 0)
            if ff & gg & 1:
                continue
            other = (nn//p, aa//p, bb//p, ff//p, gg//p)
            if whole:
                reverse.setdefault(other, []).append((state, digit))
            if other not in seen:
                seen.add(other)
                pending.append(other)
    if whole:
        return start, reverse, terminal
    return False


def binary_residues(k, c, q, modulus=45):
    # Read all accepted words backwards. The suffix grows by y -> 2y + digit.
    starts, reverse, terminals = graph(k, c, q, 2, whole=True)
    pending = deque((state, 0) for state in terminals)
    seen = set(pending)
    P = 1 << k
    while pending:
        state, y = pending.popleft()
        if state in starts:
            v = (q + 2*P*y) % modulus
            if P*v % 15 == 7 and (5*P*v+1) % 9:
                return True
        for other, digit in reverse.get(state, []):
            pair = other, (2*y+digit) % modulus
            if pair not in seen:
                seen.add(pair)
                pending.append(pair)
    return False


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("last", type=int, nargs="?", default=35)
    args = parser.parse_args()
    if not 3 <= args.last <= 35:
        parser.error("the checked table covers exponents 3 through 35")
    survivors = []
    total = 0
    for k in range(3, args.last+1):
        remaining = templates(k)
        total += len(remaining)
        output = subprocess.check_output([args.executable.resolve(), str(k), str(k), "--witnesses"], text=True)
        for line in output.splitlines():
            fields = line.split()
            if fields[0] == "excluded":
                _, exponent, c, p = fields
                assert int(exponent) == k and int(p) in (2, 3, 5)
                c, p = int(c), int(p)
                assert not graph(k, c, remaining.pop(c), p)
            elif fields[0] == "unresolved":
                _, exponent, c, q = fields
                assert int(exponent) == k
                c, q = int(c), int(q)
                assert remaining.pop(c) == q
                if binary_residues(k, c, q):
                    survivors.append((k, c, q))
        assert not remaining
        print("verified exponent", k, flush=True)
    expected = [(31, 4192014075, 741526579), (32, 6285315, 5727322795),
                (34, 34312377675, 429493347)]
    assert survivors == [row for row in expected if row[0] <= args.last]
    print(total, "templates;", "remaining", survivors)


if __name__ == "__main__":
    main()
