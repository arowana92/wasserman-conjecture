"""Compare two full cofactor enumerations and check every surviving exponent class."""
from argparse import ArgumentParser
from pathlib import Path
from subprocess import check_output

from check_cross_parts import divisors
from search_even_endpoint import prime_powers, prime_sieve


def read_report(executable, bits, marker):
    lines=check_output([str(executable.resolve()),str(bits)],text=True).splitlines()
    records=[tuple(map(int,s.split()[1:])) for s in lines if s.startswith(marker+' ')]
    assert len(records)==len(set(records))
    return set(records),tuple(map(int,lines[-1].split()))


def full_candidates(d,h):
    q,rem=divmod((1<<h)+1,d)
    if rem or q<=5 or q%8!=3 or (3*q)&(2*q) or bin(q).count('1')<7:
        return None
    first=next((i for i in range(4) if (q<<i)%15==7),None)
    if first is None or len(prime_powers(q))!=1:
        return None
    classes=[]
    for k in range(first,4*h,4):
        residue=pow(2,-k,q)
        z=residue+q if residue%2 else residue
        if z%15==0 and not z & ~(2*q):
            classes.append((k,z))
    return q,classes


def audit(bits, first_program, second_program):
    hits,counts=read_report(first_program,bits,'hit')
    long,other=read_report(second_program,bits,'long')
    assert counts[:2]==other[:2] and counts[2]==len(hits) and other[2]==len(long)
    bound=1<<bits
    recovered=set()
    small=[]
    for h in range(1,128,2):
        total=(1<<h)+1
        for d in divisors(total):
            if d==3 or d%24!=3 or d>=bound:
                continue
            q=total//d
            if h>=bits and not (3*q)&(2*q):
                recovered.add((d,h))
            if h<bits:
                small.append((d,h))

    last=127
    for (d,) in long:
        # Inverting at a wider modulus gives the first overlapping bit directly.
        size=256
        while True:
            mask=(1<<size)-1
            q0=pow(d,-1,1<<size)
            conflict=(3*q0)&(2*q0)&mask
            if conflict:
                break
            size*=2
        failure=(conflict & -conflict).bit_length()-1
        assert failure>=128
        last=max(last,failure)
        for h in range(129,failure+1,2):
            q,rem=divmod((1<<h)+1,d)
            if rem==0 and not (3*q)&(2*q):
                recovered.add((d,h))
    assert recovered==hits,(recovered-hits,hits-recovered)

    survivors=[]
    for d,h in sorted(set(small)|hits):
        result=full_candidates(d,h)
        if result is not None:
            q,classes=result
            survivors.append((d,h,q))
            if classes:
                raise ArithmeticError(('candidate exponent classes',d,h,q,classes))

    # These prime powers survive the preliminary tests and occur at nonminimal h.
    print('Distinct admissible Q before the inverse test:',sorted({q for _,_,q in survivors}))
    print('Prefix counts:',counts,'independent long-prefix counts:',other)
    print('Latest independently detected overlap for long prefixes:',last)
    print('All modular exponent classes excluded for',len(survivors),'remaining (d,h,Q) triples.')
    if bits==52:
        assert counts==(879862521,293293896,6,186)
        assert other==(879862521,293293896,18)


def radical_check():
    exceptional=[p for p in prime_sieve(400000) if pow(2,p-1,p*p)==1]
    assert exceptional==[1093,3511]
    for p,order in ((1093,364),(3511,1755)):
        assert pow(2,order,p)==1
        assert all(pow(2,order//r,p)!=1 for r in prime_powers(order))
        assert order%4!=2
    assert (164*400000)**2 < 2**52
    print('Checked all odd primes through 400,000 for the radical bound.')


if __name__=='__main__':
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('primary',type=Path)
    parser.add_argument('cross',type=Path)
    parser.add_argument('--bits',type=int,default=52)
    args=parser.parse_args()
    radical_check()
    audit(args.bits,args.primary,args.cross)
