"""Finite direct orbits for the independent coordinates of the cube-base sieve."""
from fractions import Fraction
from math import gcd
from check_cofactor_986_powers_cross import extendible_prefix


def logs_29():
    modulus = 29**3
    x = 1
    logs = {}
    while x not in logs:
        logs[x] = len(logs)
        x = x*2%modulus
    assert len(logs) == 28*29**2
    return logs


def binary_orbit(q,n):
    w = pow(q,n,1024)
    step = pow(q,84,1024)
    start = w
    while True:
        v = (986*w+1)*pow(7,-1,2048)%2048
        if (v+(2*v%2048)+(4*v%2048) == (v | (2*v%2048) | (4*v%2048))):
            return True
        w = w*step%1024
        if w == start:
            return False


def verify():
    logs = logs_29()
    period = len(logs)
    total = Fraction()
    for slopes,target in (((119,203,664),10),((203,714,69),21)):
        eligible = [x for x in range(target,29**3,29)
                    if extendible_prefix(slopes,x,29,3)]
        targets = [logs[x] for x in eligible]
        mass = Fraction()
        for base in range(1,29):
            if len({pow(base,k,29) for k in range(28)}) != 28:
                continue
            n, = [k for k in range(3,85,3) if pow(base,k,29)==target]
            good29 = 0
            for q in range(base,29**3,29):
                i = logs[q]
                step = gcd(84*i,period)
                good29 += any((j-n*i)%step==0 for j in targets)
            good7 = sum(q%7!=0 and extendible_prefix(slopes,pow(q,n,49),7,2)
                        and pow(q,n,7)==1 for q in range(1,49))
            good17 = sum(any(pow(q,n+84*k,289)==84 for k in range(17))
                         for q in range(1,289) if q%17)
            good2 = sum(binary_orbit(q,n) for q in range(1,1024,2))
            good15 = sum(gcd(q,15)==1 and gcd(986*pow(q,n,15)+1,15)==1 for q in range(15))
            mass += (Fraction(good2,512)*Fraction(good7,42)*Fraction(good17,272)
                     *Fraction(good29,period)*Fraction(good15,8))
        print(slopes,"orbit count",mass,flush=True)
        total += mass
    assert total == Fraction(99,188384)
    print("Independent CRT density",total)


if __name__ == '__main__':
    verify()
