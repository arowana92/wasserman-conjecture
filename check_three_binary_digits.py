"""Check the short-relation bounds and residue divisibilities on sample exponent pairs."""

from math import gcd, isqrt


def short_relation(s, t):
    limit = isqrt(t - 1) + 1
    for b in range(1, limit + 1):
        a = (2*b*s + t) // (2*t)
        if 1 <= a <= b and abs(b*s - a*t)**2 <= t:
            d = gcd(a, b)
            return a // d, b // d
    raise ArithmeticError(("missing approximation", s, t))


def check(s, t):
    h_limit = 300
    if not (s*s > t and t-s >= 4):
        raise ValueError("pair lies outside the short-relation step")
    a, b = short_relation(s, t)
    delta = b*s - a*t
    left = [(u**b) << max(-delta, 0) for u in range(h_limit)]
    right = [(v**a) << max(delta, 0) for v in range(h_limit)]
    x, y = 2**s, 2**t
    omitted = set()
    checked = 0
    for u in range(h_limit - 1):
        for v in range(h_limit - 1 - u):
            relation = left[u] - right[v]
            modulus = gcd(x - u, y - v)
            if modulus and relation % modulus:
                raise ArithmeticError((s, t, u, v, modulus))
            if relation == 0:
                omitted.add(1 + u + v)
            checked += 1
    degree = sum(h_limit-h for h in range(h_limit) if h not in omitted)
    if len(omitted) > 18 or degree < 39750:
        raise ArithmeticError((s, t, omitted, degree))
    return checked, len(omitted), degree


if __name__ == "__main__":
    pairs = [(12, 16), (15, 30), (22, 43), (19, 39), (73, 119),
             (229, 541), (401, 983), (987, 1597)]
    total = 0
    for s, t in pairs:
        checks, omitted, degree = check(s, t)
        print(f"s={s}, t={t}, omitted offsets={omitted}, remaining degree={degree}")
        total += checks
    print(f"Checked {total} residue pairs.")
