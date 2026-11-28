"""Check the logarithmic block construction with free high digits at unpaired primes."""

from math import gcd, isqrt, log, prod, sqrt

from check_consecutive_certificate_lower_bound import next_power, primes_through


def cut_unit(total, part, p):
    while part:
        part,digit = divmod(part,p)
        total,available = divmod(total,p)
        if digit > available:
            return False
    return True


def factorial_floor_unit(n, h, a, p):
    large = n-h-a
    while h or a or n != large:
        n,h,a,large = n//p,h//p,a//p,large//p
        if n != h+a+large:
            return False
    return True


def block(length, cutoff, exponent):
    t = length+1
    assert 3 <= t <= cutoff
    primes = primes_through(cutoff)
    small = max(isqrt(t)+1,int(t/log(t)**4))
    middle = max(small,int(t/sqrt(log(t))))
    q = {p:next_power(p,t) for p in primes if p <= t}
    needed = dict(q)
    groups = {}
    paired = set()
    for h in range(2,t+1):
        active = [p for p in q if h*(h-1) % p == 0]
        left = [p for p in active if p <= small]
        right = [p for p in active if p > small]
        assert len(right) <= 2
        if len(right) == 2:
            paired.update(right)
        for group in (left,right):
            a = prod(q[p] for p in group)
            for p in group:
                needed[p] = max(needed[p],a)
        groups[h] = left,right
    free = {p for p in q if p > middle and p not in paired}
    base,modulus = 0,1
    for p in primes:
        if p in q and p not in free:
            m, residue = next_power(p,needed[p]), -q[p]
        else:
            m, residue = p,0
        base += modulus*((residue-base)*pow(modulus,-1,m) % m)
        modulus *= m
    start = 2**exponent
    base = start+(base-start)%modulus
    tried = 1
    while any(base//p%p+t//p >= p for p in free):
        base += modulus
        tried += 1
        assert tried <= t**3
    primorial = prod(primes)
    unpaired = 0
    for h,(left,right) in groups.items():
        n = base+h
        assert n%h == (n-1)%(h-1) == 0
        assert gcd(n*(n-1)//(h*(h-1)),primorial) == 1
        for group in (left,right):
            if group and group[0] in free:
                p, = group
                a = p
                while base % (p*a) == 0:
                    a *= p
                unpaired += 1
            else:
                a = prod(q[p] for p in group)
            assert 0 < a < base
            assert all(cut_unit(n,h,p) and cut_unit(base,a,p) for p in group)
            assert all(factorial_floor_unit(n,h,a,p) for p in group)
        assert (n-start)**7 <= 6**6*(n-1)**4 and n < 2*start
    print("length",length,"cutoff",cutoff,"binary exponent",exponent,"modulus bits",
          modulus.bit_length(),"tested origins",tried,"offset bits",(base-start).bit_length())
    print("small cutoff",small,"middle cutoff",middle,"paired primes",len(paired),
          "free primes",len(free),"singleton splits",unpaired)


if __name__ == "__main__":
    block(128,1000,4096)
    block(1000,1009,16384)
