"""Check a minimum four-prime certificate using standard-library integer arithmetic."""
from math import gcd, isqrt

# For each nontrivial node, these primes supply enough of n-1 for
# Pocklington's criterion. The witness is 2 throughout the tree.
FACTORS = {
    324994591: (521, 239),
    2771149619312071357067: (324994591, 303013),
    3226141379: (22807, 661),
    3329377903129: (3226141379,),
    31005227028198764011132356842134238486295073: (2771149619312071357067, 3329377903129),
    10039633: (209159,),
    7148218697: (10039633,),
    921810807957433055393: (7148218697, 502669),
    519901295687992243241653: (921810807957433055393,),
    1814839290033733262260073929992575167147: (519901295687992243241653,),
    22832816491: (739, 103, 101),
    10325410607927039: (22832816491,),
    1621418506274976531096954667: (10325410607927039,),
    155656176602397746985307648033: (1621418506274976531096954667,),
    205407527369917799624694679204427930885246377: (155656176602397746985307648033,),
    2993155353253694493393129679484290878994324303904831: (205407527369917799624694679204427930885246377,),
}


def prove_prime(n):
    if n < 10**6:
        assert n >= 2 and all(n % j for j in range(2, isqrt(n)+1))
        return
    assert pow(2, n-1, n) == 1
    supported = 1
    for q in FACTORS[n]:
        prove_prime(q)
        assert (n-1) % q == 0
        assert gcd(pow(2, (n-1)//q, n)-1, n) == 1
        cofactor = n-1
        while cofactor % q == 0:
            cofactor //= q
            supported *= q
    assert (n-1) % supported == 0 and supported**2 > n


if __name__ == '__main__':
    heights = [13*k+4*k*k for k in range(6)]
    m = sum(1 << k for k in heights)
    n = 1+64*m
    q, s, r = n//96537121, m//(3*(2**33+1)), n-2
    assert q*96537121 == n and s*3*(2**33+1) == m
    assert 1 < n//q < q and len({2,q,s,r}) == 4
    for prime in (2,q,s,r):
        prove_prime(prime)
    assert n-(1 << (6+heights[-1])) < min(2*q, s)
    assert n < q*(2*s+1) and n-2 == r and n % 2
    print(n, (2,q,s,r))
