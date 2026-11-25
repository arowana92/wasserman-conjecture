"""Finite carry graphs for rows with entries affine in a nonnegative input."""

class CarryGraph:
    def __init__(self, base, slopes, offsets):
        self.base = base
        self.slopes = slopes
        self.offsets = offsets
        self.states = [tuple(offsets)]
        index = {self.states[0]: 0}
        self.edges = []
        for state in self.states:
            successors = []
            for digit in range(base):
                values = [s + a*digit for s, a in zip(state, slopes)]
                if sum(x % base for x in values) >= base:
                    continue
                destination = tuple(x//base for x in values)
                if destination not in index:
                    index[destination] = len(self.states)
                    self.states.append(destination)
                successors.append((digit, index[destination]))
            self.edges.append(successors)
        self.accept = sum(1 << i for i, row in enumerate(self.states) if self.no_carry(row))

    def no_carry(self, entries):
        row = tuple(entries)
        while any(row):
            if sum(x % self.base for x in row) >= self.base:
                return False
            row = tuple(x//self.base for x in row)
        return True

    def block_bound(self, length):
        # One path per starting state and readable word, allowing every endpoint.
        ways = [1]*len(self.states)
        for _ in range(length):
            ways = [sum(ways[j] for _, j in successors) for successors in self.edges]
        return sum(ways)

    def accepts_integer(self, t):
        return self.no_carry([a*t+b for a, b in zip(self.slopes, self.offsets)])


def forced_bases(shape, primes_in_d, primes_in_u):
    # Avoidance in C at p | d requires two zero residues in its normalized row.
    ys = shape[4:]
    return sorted({p for p, _ in primes_in_u} |
                  {p for p, w in primes_in_d if sum(y % w == 0 for y in ys) < 2})


class PrefixMatcher:
    """Read high digits backward through a low-digit-first carry graph."""

    def __init__(self, graph, limit):
        self.graph = graph
        self.reverse = [{} for _ in range(graph.base)]
        for source, edges in enumerate(graph.edges):
            for digit, target in edges:
                by_target = self.reverse[digit]
                by_target[target] = by_target.get(target, 0) | (1 << source)
        self.cache = {}
        self.powers = [1]
        self.reachable = [1]  # Only the initial state before reading any low digits.
        fixed = False
        while self.powers[-1] <= limit:
            self.powers.append(self.powers[-1]*graph.base)
            mask = self.reachable[-1]
            if fixed:
                self.reachable.append(mask)
                continue
            following = 0
            for source, edges in enumerate(graph.edges):
                if mask & (1 << source):
                    for _, target in edges:
                        following |= 1 << target
            self.reachable.append(following)
            fixed = following == mask

    def previous(self, mask, digit):
        key = (mask, digit)
        if key not in self.cache:
            table = self.reverse[digit]
            previous = 0
            while mask:
                low = mask & -mask
                previous |= table.get(low.bit_length()-1, 0)
                mask -= low
            self.cache[key] = previous
        return self.cache[key]


def intersection(graphs, limit, modulus=1, residue=0):
    """Return all inputs in [2,limit] in the given congruence class passing every graph."""
    matchers = [PrefixMatcher(g, limit) for g in graphs]
    # Each cell represents [origin,origin+base**remaining). The mask contains
    # states from which its fixed high digits can be read to acceptance.
    cells = tuple((len(m.powers)-1, 0, m.graph.accept) for m in matchers)
    stack = [(2, limit, cells)]
    survivors = []
    nodes = tested = 0
    while stack:
        lo, hi, cells = stack.pop()
        nodes += 1
        if hi-lo < 32:
            for t in range(lo+(residue-lo) % modulus, hi+1, modulus):
                tested += 1
                if all(g.accepts_integer(t) for g in graphs):
                    survivors.append(t)
            continue

        i = max(range(len(cells)), key=lambda j: matchers[j].powers[cells[j][0]])
        m = matchers[i]
        remaining, origin, mask = cells[i]
        width = m.powers[remaining-1]
        for digit in range((lo-origin)//width, (hi-origin)//width+1):
            active = m.previous(mask, digit) & m.reachable[remaining-1]
            if not active:
                continue
            start = origin+digit*width
            replacement = (remaining-1, start, active)
            children = cells[:i]+(replacement,)+cells[i+1:]
            stack.append((max(lo, start), min(hi, start+width-1), children))
    return sorted(survivors), nodes, tested
