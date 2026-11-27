"""Exact congruence feasibility for a sequence of candidate exponent classes."""
from itertools import combinations
from math import gcd
from typing import NamedTuple


class Candidate(NamedTuple):
    d: int
    a: int
    h: int
    forbidden_length: int
    ternary: bool = False


SMALL = (Candidate(12031,769251571,2043279108,3,True),
         Candidate(56279,801,1743,4),Candidate(56279,717,1743,4),
         Candidate(60829,35049,627270,4),Candidate(61423,7030,9520,4),
         Candidate(65095,40240,80500,3),Candidate(65243,174688,287820,3,True))
LARGE = (Candidate(7338895,2337019,2935556,3,True),
         Candidate(14117615,3948614,5647044,3,True))


def meeting(first, stride, other, period):
    # A zero stride specifies a single integer.
    if stride == 0:
        return (first,0) if (first-other)%period == 0 else None
    g = gcd(stride,period)
    if (other-first)%g:
        return None
    q = period//g
    index = ((other-first)//g * pow(stride//g,-1,q)) % q
    return first+stride*index,stride*q


def allowed_indices(indices, row):
    last = indices[-1]
    for step in range(1,last//(row.forbidden_length-1)+1):
        if all(last-i*step in indices for i in range(1,row.forbidden_length)):
            return False
    if row.ternary:
        for i,j in combinations(indices[:-1],2):
            x,y = last-i,last-j
            while x%3 == y%3 == 0:
                x //= 3
                y //= 3
            if x%3 and y%3 and (x-y)%3:
                return False
    return True


def extend(prefix, position, row, step_class):
    earlier = [i for i,r in enumerate(prefix) if r == row]
    if not allowed_indices(earlier+[position],row):
        return None
    for i,old in enumerate(prefix):
        m = gcd(row.h,old.h)
        diff = row.a-old.a
        gap = position-i
        g = gcd(gap,m)
        if diff%g:
            return None
        modulus = m//g
        residue = (diff//g*pow(gap//g,-1,modulus)) % modulus
        step_class = meeting(*step_class,residue,modulus)
        if step_class is None:
            return None
    return step_class


def lengths(rows, limit, step=None):
    active = [((),(0,1) if step is None else (step,0))]
    counts = []
    for position in range(limit):
        after = []
        for prefix,step_class in active:
            for row in rows:
                new_step = extend(prefix,position,row,step_class)
                if new_step is not None:
                    after.append((prefix+(row,),new_step))
        if rows == SMALL and step is None and position == 3:
            assert {c for _,c in after} == {(h,420) for h in (125,155,265,295)}
        counts.append(len(after))
        active = after
    assert not active
    return counts


if __name__ == "__main__":
    for rows,limit,step in ((SMALL,5,None),(SMALL,3,1),(LARGE,3,None),
                            (SMALL+LARGE,12,None),(SMALL+LARGE,4,1)):
        print(len(rows),step,lengths(rows,limit,step))
