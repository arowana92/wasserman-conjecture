"""Independent small-product audit by factoring m^2+L. Requires python-flint."""
import argparse
from math import gcd, lcm
from check_cross_cut_audit import common_primes, divisors


def primitive_product(n, a):
    d = gcd(a, n)
    return (a//d)*((n-a)//d)


def scan(limit, bound):
    lowest = parity = 0
    survivors, critical = [], []
    for ell in range(1, bound+1):
        for m in range(1, limit):
            for x in divisors(m*m+ell):
                if not m < x <= limit:
                    continue
                for h in range(1, bound//ell+1):
                    y, rem = divmod(m*m*(x-m)-h, ell)
                    if rem or y < x or gcd(x,y) != 1:
                        continue
                    pred, rem = divmod(m*(x-m)*(y-m), h)
                    if rem:
                        continue
                    n = pred+1
                    g, rem = divmod(m*n, x*y)
                    if rem or g < 1:
                        continue
                    a, b = g*x, g*y
                    if lcm(primitive_product(n,a), primitive_product(n,b)) % pred:
                        continue
                    lowest += 1
                    if a & (n-a) and b & (n-b):
                        continue
                    parity += 1
                    if not any(common_primes(n, a, b)):
                        survivors.append((h,ell,x,m,n,a,b))
                    if h == ell == 1:
                        layer = 3
                        while layer <= n and min(a%layer, b%layer) <= n%layer:
                            layer *= 3
                        critical.append((x,m,n,layer if layer <= n and n%3 <= 1 else None))
    return lowest, parity, sorted(survivors), sorted(critical)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("limit", type=int, nargs="?", default=10**6)
    parser.add_argument("product", type=int, nargs="?", default=1)
    options = parser.parse_args()
    result = scan(options.limit, options.product)
    print("Lowest layers, parity, endpoint-coprime pairs:", *result[:3])
    print("h=L=1 parity survivors with base-three carry layers:", result[3])
    if (options.limit, options.product) == (10**6, 1):
        assert result == (121202, 2, [],
                          [(22405, 17648, 124380638800347396945, 9),
                           (24565, 2928, 11751878708987668945, 27)])
