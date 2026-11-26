"""Independent integer-language and modular checks at denominators 364 and 880."""
from collections import deque
from math import gcd

from carry_blocks import CarryGraph, accepted_hull


def binary_certificate(slopes, power, low, modulus):
    entries = [a*low+b for a,b in zip(slopes, (0,0,1))]
    assert sum(x % power for x in entries) == (entries[0] | entries[1] | entries[2]) % power
    graph = CarryGraph(2, slopes, [x//power for x in entries])
    reference, period = accepted_hull(graph)
    first, step = low+power*reference, power*period
    assert CarryGraph(2, slopes, (0,0,1)).accepts_integer(first)

    # Forward search, retaining the absolute integer residue and place. Test
    # the raw column bits and terminal integers directly, without CarryGraph.
    initial = tuple(x//power for x in entries), low % modulus, power % modulus
    seen = {initial}
    todo = deque([initial])
    residues = set()
    while todo:
        state, value, place = todo.popleft()
        if sum(state) == (state[0] | state[1] | state[2]):
            residues.add(value)
        for bit in (0,1):
            x,y,z = (c+a*bit for a,c in zip(slopes,state))
            if ((x&y) | (x&z) | (y&z)) & 1:
                continue
            nxt = (x//2,y//2,z//2), (value+bit*place)%modulus, 2*place%modulus
            if nxt not in seen:
                seen.add(nxt)
                todo.append(nxt)
    print(slopes, "hull", first % step, "mod", step, "forward residues", sorted(residues))
    assert residues == {first % modulus} and step % modulus == 0
    return first % step, step


if __name__ == '__main__':
    assert binary_certificate((76,247,41),4,3,135) == (135,69120)
    for slopes in ((76,247,41),(5,7,7)):
        for digit in (1,2):
            assert sum(a*digit%3 for a in slopes) >= 3
    assert binary_certificate((272,561,47),16,1,5) == (241,15360)
    for w in range(1,17*25,5):
        if (880*w+1)%17:
            continue
        v = (880*w+1)//17
        assert v%25 == 18 and sum(a*v%25 for a in (2,5,10)) == 31
    # A nontrivial gcd of input differences and a singleton are distinguished.
    assert accepted_hull(CarryGraph(2,(1,1),(0,0))) == (0,0)
    assert gcd(*accepted_hull(CarryGraph(2,(1,0),(0,0)))) == 1
