"""Verify exponent progression bounds by exhausting two small residue rings."""
from functools import lru_cache
from itertools import combinations


# (cofactor, necessary residue, modulus, excluded progression length, ternary test)
SMALL = ((12031,7,12,3,True),(56279,3,21,4,False),(60829,189,210,4,False),
         (61423,30,140,4,False),(65095,60,140,3,False),(65243,28,60,3,True))
LARGE = ((7338895,3,4,3,True),(14117615,2,12,3,True))


def equal_ternary_valuations(i,j,k):
    a,b = j-i,k-j
    while a%3 == b%3 == 0:
        a //= 3
        b //= 3
    return a%3 == b%3 != 0


def steps(rows, period, length, step_one=False):
    choices = []
    from math import gcd
    for e in range(period):
        choices.append(tuple(i for i,(_,a,h,_,_) in enumerate(rows)
                             if (e-a)%gcd(h,period) == 0))
    if rows == SMALL:
        assert all(len(c) <= 1 for c in choices)
        assert sum(bool(c) for c in choices) == 70

    forbidden = []
    for _,_,_,size,ternary in rows:
        sets = []
        for stride in range(1,length):
            for i in range(length-(size-1)*stride):
                sets.append(sum(1 << (i+j*stride) for j in range(size)))
        if ternary:
            sets.extend((1<<i)+(1<<j)+(1<<k) for i,j,k in combinations(range(length),3)
                        if equal_ternary_valuations(i,j,k))
        forbidden.append(sets)

    @lru_cache(None)
    def possible(word):
        marks = [0]*len(rows)
        def assign(i):
            if i == length:
                return True
            for r in word[i]:
                marks[r] |= 1 << i
                if all(marks[r]&m != m for m in forbidden[r]) and assign(i+1):
                    return True
                marks[r] ^= 1 << i
            return False
        return assign(0)

    surviving = set()
    for h in ((1,) if step_one else range(period)):
        for b,first in enumerate(choices):
            if not first:
                continue
            word = tuple(choices[(b+j*h)%period] for j in range(length))
            if all(word) and possible(word):
                surviving.add(h)
                break
    return surviving


if __name__ == "__main__":
    assert steps(SMALL,420,4) == {19,41,125,155,265,295,379,401}
    assert not steps(SMALL,420,5) and not steps(SMALL,420,3,True)
    assert not steps(SMALL+LARGE,84,12) and not steps(SMALL+LARGE,84,4,True)
    print("All starting residues and steps verified modulo 420 and 84.")
