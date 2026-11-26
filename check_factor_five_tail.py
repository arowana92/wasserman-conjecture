"""Certify the two exceptional exponent progressions by ternary paths and quinary carries."""

from collections import Counter, defaultdict, deque
from typing import NamedTuple

from flint import fmpz

from check_factor_five_rays import initial_states, no_carries, successors


class Progression(NamedTuple):
    q: int
    z: int
    first: int
    stride: int

    @property
    def slopes(self):
        return 3*self.q, self.z, 2*self.q-self.z


RAYS = (Progression(1467779,2918400,2337019,2935556),
        Progression(2823523,5636100,3948614,5647044))
DEPTH = 256
BLOCK = 2048

# Reachable, terminal, productive, productive outdegree histogram, maximum block count.
GRAPH_COUNTS = {
    1467779: (151397,21,29656,{1:28350,2:1306},40879096506588498826141913655374463598481),
    2823523: (199807,46,36778,{1:35089,2:1689},878909565367958979420167123851270494639),
}
BASE_FIVE_COUNTS = {
    1467779: {2:33467,3:46664,4:2073,5:1035,6:213,7:13,8:35,10:2},
    2823523: {2:81956,3:82110,4:26400,5:9162,6:4445,7:881,8:160,9:34,10:13},
}


def graph_bounds(ray, valuation):
    initial = list(initial_states(ray.q,ray.z,3,ray.first,ray.stride))
    slopes = ray.slopes
    queue = deque(initial)
    visited = set(initial)
    terminal = []
    edges = {}
    parents = defaultdict(list)
    while queue:
        node = queue.popleft()
        if no_carries(node,3):
            terminal.append(node)
        edges[node] = list(successors(node,slopes,3))
        for child in edges[node]:
            parents[child].append(node)
            if child not in visited:
                visited.add(child)
                queue.append(child)

    # A finite integer input can use only states from which termination is possible.
    productive = set(terminal)
    queue = deque(terminal)
    while queue:
        for parent in parents[queue.popleft()]:
            if parent not in productive:
                productive.add(parent)
                queue.append(parent)
    index = {v:i for i,v in enumerate(productive)}
    targets = [[index[w] for w in edges[v] if w in productive] for v in productive]
    degrees = Counter(map(len,targets))
    assert set(degrees) == {1,2}
    first = [row[0] for row in targets]
    second = [(i,row[1]) for i,row in enumerate(targets) if len(row)==2]
    paths = [1]*len(targets)
    for _ in range(BLOCK):
        after = [paths[j] for j in first]
        for i,j in second:
            after[i] += paths[j]
        paths = after
    maximum = max(paths)
    assert (len(visited),len(terminal),len(productive),degrees,maximum) == GRAPH_COUNTS[ray.q]
    assert maximum**625 < 3**(26*BLOCK)
    if ray.q == 2823523:
        assert maximum**25 < 3**BLOCK
    print('Q',ray.q,'maximum 2048-digit path count',maximum)

    # Count all low words independently, including states unable to terminate.
    paths = Counter(initial)
    initial_digits = max(valuation,2)
    for _ in range(DEPTH+valuation-initial_digits):
        after = Counter()
        for node,count in paths.items():
            for w in edges[node]:
                after[w] += count
        paths = after
    return sum(paths.values())


def finite_range(ray, valuation):
    fixed_modulus = 3**valuation
    modulus = fixed_modulus*3**DEPTH
    first = pow(2,ray.first,modulus)
    stride = pow(2,ray.stride,modulus)
    a,b,c = ray.slopes
    assert no_carries(((a*first)%fixed_modulus,(b*first)%fixed_modulus,
                       (c*first+1)%fixed_modulus),3)
    states = [(0,first)]
    place = 1
    for j in range(DEPTH):
        new = []
        layer = 3*fixed_modulus*place
        square = stride*stride % modulus
        for t,w in states:
            for digit,x in enumerate((w,w*stride%modulus,w*square%modulus)):
                if j==0 and (5*ray.q*x+1)%9==0:
                    continue
                if (a*x)%layer + (b*x)%layer + (c*x+1)%layer < layer:
                    new.append((t+digit*place,x))
        states = new
        place *= 3
        stride = pow(stride,3,modulus)

    assert len(states) == len({t for t,_ in states})
    rejected = Counter()
    for t,w in states:
        assert 0<=t<3**DEPTH
        k = ray.first+ray.stride*t
        # Recalculate the residues from the exponent via an independent library.
        assert pow(fmpz(2),k,modulus) == w
        assert no_carries(((a*w)%modulus,(b*w)%modulus,(c*w+1)%modulus),3)
        v = pow(2,k,5**10)
        for j in range(1,11):
            size = 5**j
            if (5*ray.q*v+1)//size > sum(x//size for x in (a*v,b*v,c*v+1)):
                rejected[j] += 1
                break
        else:
            raise ArithmeticError(('requires more digits',ray.q,t))
    assert rejected == BASE_FIVE_COUNTS[ray.q]
    assert ray.first+ray.stride*3**DEPTH > 10**128
    print('Q',ray.q,'first overflowing base-five layers',dict(sorted(rejected.items())))
    return len(states)


if __name__=='__main__':
    for ray in RAYS:
        valuation = 1
        while pow(2,ray.stride,3**(valuation+1))==1:
            valuation += 1
        assert pow(2,ray.stride,3**valuation)==1
        assert graph_bounds(ray,valuation) == finite_range(ray,valuation)
