"""Check the progression exclusions using rational cells and direct scalar residues."""
from math import lcm

from check_binary_exponent_progressions import published_rays


def cell_graph(slopes, p, first_digit):
    # The three floors at a common real y change only at these breakpoints.
    # A shift of one in the last numerator can add at most one to its floor.
    scale = lcm(*slopes)
    cuts = {i*(scale//s) for s in slopes for i in range(s)}
    states = set()
    x, y, z = slopes
    for cut in cuts:
        a, b, c = (s*cut//scale for s in slopes)
        states.add((a, b, c))
        states.add((a, b, c+1))
    states = sorted(states)
    index = {s: i for i, s in enumerate(states)}
    parents = [[] for _ in states]
    edges = [[] for _ in states]
    for i, (a, b, c) in enumerate(states):
        for digit in range(p):
            A, B, C = a+x*digit, b+y*digit, c+z*digit
            target = A//p, B//p, C//p
            if sum(target) == (A+B+C)//p:
                child = index[target]
                parents[child].append(i)
                edges[i].append((digit, child))
    zero = index[0, 0, 0]
    complete, work = {zero}, [zero]
    for j in work:
        for i in parents[j]:
            if i not in complete:
                complete.add(i)
                work.append(i)
    root = dict(edges[index[0, 0, 1]])[first_digit]
    assert root in complete
    reached, work = {root}, [root]
    output = {}
    for i in work:
        output[i] = {digit: j for digit, j in edges[i] if j in complete}
        for j in output[i].values():
            if j not in reached:
                reached.add(j)
                work.append(j)
    return output, root


def check_extensions(edges, p, terms, state, depth, prefix=0, prefix_power=1):
    # In the linear test lift a starting block and its unit difference.
    # In the geometric test lift the first scalar and the ratio 1+prefix_power*b.
    def visit(k, width, a, b, states):
        if k == depth:
            raise ArithmeticError(("progression survived", prefix, a, b))
        following = p*width
        for digit, child in edges[states[0]].items():
            first = a+width*digit
            for step_digit in range(p):
                if k == 0 and step_digit == 0:
                    continue
                step = b+width*step_digit
                word = first
                next_states = [child]
                for s in states[1:]:
                    if prefix_power == 1:
                        word = (word+step) % following
                    else:
                        # Divide out the common low prefix in q -> q*(1+p^t*step).
                        word = (word + (prefix+prefix_power*word)*step) % following
                    dest = edges[s].get(word//width)
                    if dest is None:
                        break
                    next_states.append(dest)
                else:
                    visit(k+1, following, first, step, next_states)

    visit(0, 1, 0, 0, [state]*terms)


# The linear depths and uniform maximum geometric depths reported by the first verifier.
DEPTHS = ((1,0), (12,9), (1,0), (13,5), (8,5), (5,1), (8,4),
          (5,1), (4,1), (12,4), (1,0))


def main():
    for ray, (depth, shallow) in zip(published_rays(), DEPTHS):
        row, p, terms = ray
        first = pow(2, row.a, p)
        assert pow(2, row.h, p) == 1
        edges, root = cell_graph(row.slopes, p, first)
        for state in edges:
            check_extensions(edges, p, terms, state, depth)
        prefixes = [(first, root)]
        for t in range(1, depth):
            width = p**t
            for a, state in prefixes:
                check_extensions(edges, p, terms, state, shallow, a, width)
            prefixes = [(a+digit*width, j) for a, i in prefixes for digit,j in edges[i].items()]
        print(row.d, row.a, len(edges), "cell states verified", flush=True)


if __name__ == "__main__":
    main()
