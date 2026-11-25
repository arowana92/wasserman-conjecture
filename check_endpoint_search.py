"""Compare the CRT row search against exact binomial products at small weights."""

from math import comb
from search_endpoint import atoms, find_row, smallest_factors, splits


def main():
    spf = smallest_factors(120)
    tested = feasible = 0
    for n in range(3, 121):
        values = [comb(n, a)*comb(n-a, b) for a in range(1, n//3+1)
                  for b in range(a, (n-a)//2+1)]
        for u, u_factors, _ in splits(atoms(n, spf)):
            for s, s_factors, _ in splits(atoms(n-1, spf)):
                primes = [p for p, _ in u_factors+s_factors]
                possible = any(all(value % p for p in primes) for value in values)
                row = find_row(n, u, u_factors, s, s_factors)
                if possible != bool(row):
                    raise ArithmeticError((n, u, s, possible, row))
                if row:
                    a, b, c = row
                    value = comb(n, a)*comb(n-a, b)
                    if min(row) < 1 or a+b+c != n or any(value % p == 0 for p in primes):
                        raise ArithmeticError(row)
                tested += 1
                feasible += possible
    print(f"Matched exact enumeration on {tested} assignments, {feasible} feasible.")


if __name__ == "__main__":
    main()
