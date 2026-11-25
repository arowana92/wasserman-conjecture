"""Check the large prime-power layers of the weighted lemma with exact coefficients."""

from itertools import permutations
from math import comb, isqrt, prod


def factor_product(x, t):
    return prod((x-j)**(t-j) for j in range(t))


def main():
    checks = nontrivial = 0
    primes = [p for p in range(2, 81) if all(p % d for d in range(2, isqrt(p)+1))]
    orders = sorted(set(order for shape in [(2, 2, 2), (1, 5, 0), (2, 4, 0),
                                             (1, 2, 3), (4, 1, 1), (6, 6, 6)]
                        for order in permutations(shape)))
    for n in range(12, 81):
        for cutoff in orders:
            h_limit = sum(cutoff)
            if h_limit > n:
                continue
            # Extract precisely the powers expected from layers exceeding H.
            required = {}
            for p in primes:
                exponent = 0
                power = p
                while power <= n:
                    residue = n % power
                    if power > h_limit and residue < h_limit:
                        exponent += h_limit-residue
                    power *= p
                if exponent:
                    required[p] = p**exponent
            for a in range(1, n//3+1):
                for b in range(a, (n-a)//2+1):
                    row = a, b, n-a-b
                    if any(x < t for x, t in zip(row, cutoff)):
                        continue
                    coefficient = comb(n, a)*comb(n-a, b)
                    product = prod(factor_product(x, t) for x, t in zip(row, cutoff))
                    for p, divisor in required.items():
                        if coefficient % p:
                            checks += 1
                            nontrivial += divisor > p
                            if product % divisor:
                                raise ArithmeticError((n, row, cutoff, p))
    print(f"Passed {checks} exact prime/row/cutoff checks, {nontrivial} with divisor > p.")


if __name__ == "__main__":
    main()
