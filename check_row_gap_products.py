"""Verify the unrestricted gap identities with exact trinomial coefficients."""

from itertools import combinations, permutations
from math import gcd, prod

from check_endpoint_pair_bound import rows


def check_coefficients(limit):
    single = double = cancellation = forced = 0
    for weight in range(3,limit+1):
        table = [(tuple(b-a for a,b in combinations(row,2) if b>a),coefficient)
                 for row,coefficient in rows(weight)]
        for gaps,first in table:
            for other,second in table:
                common = gcd(first,second)
                D,J = gcd(weight,common),gcd(weight-1,common)
                for d in gaps:
                    assert D*d >= 3
                    if d >= 2:
                        scale = gcd(weight//D,d)
                        assert 27*(weight-1)**2 <= J**2*D**3*d**2*scale*(d*d-1)**2
                        single += 1
                    for e in other:
                        tau,rem = divmod(D*d*e,weight)
                        assert rem == 0
                        diff = (D*(d-e))**2-(tau-D)**2
                        assert J*tau*diff % (weight-1) == 0
                        assert J*tau*(d*d-1)*((D*d)**2-tau**2) % (weight-1) == 0
                        double += 1
                        cancellation += diff == 0
                        # Test the real-valued gcd lower bound by clearing N^3.
                        H = d*e*abs((weight*(d-e))**2-(d*e-weight)**2)
                        if H:
                            assert (D*J)**3*H >= (weight-1)*weight**3
                            forced += H < (weight-1)*weight**3
    print('Exact coefficient checks:',single,'single gaps,',double,'gap pairs;',
          cancellation,'vanishing expressions,',forced,'nonvanishing bounds above one')


def check_binary_cutoff(limit):
    checked = even_multiples = 0
    # Distribute each occupied bit of the weight to exactly one part.
    for weight in range(3,limit+1):
        parts = [(0,0,0)]
        for index in range(weight.bit_length()):
            bit = 1 << index
            if weight & bit:
                parts = [(a+bit,b,c) for a,b,c in parts] + [(a,b+bit,c) for a,b,c in parts] + [
                    (a,b,c+bit) for a,b,c in parts]
        run = max(map(len,bin(weight)[2:].split('0')))
        for a,large,third in parts:
            if not 0 < a < large or not third:
                continue
            gap = large-a
            assert a < (1 << (run+1))*gap
            denominator = gcd(a,third)
            b = (2*a+third)//denominator
            assert a < b*gap
            denominator = gcd(large,third)
            c = (2*large+third)//denominator
            if third % (2*large):
                assert large < (c+1)*gap
            else:
                even_multiples += 1
            checked += 1
    print('Binary cutoff checks:',checked,'including',even_multiples,'even multiples')


def check_low_residues():
    from check_six_part_product import assignments

    systems = set()
    for order in set(permutations((1,1,2))):
        systems.update(assignments((1,1,1),order))
    total = resonance = 0
    for system in systems:
        parts = []
        for row in system:
            modulus = prod(w for _,w in row)
            parts.append([f*(modulus//w) for f,w in row])
        a,c = parts
        U,V = sum(c),sum(a)
        N = U*V
        for i,j in combinations(a,2):
            for l,m in combinations(c,2):
                d,e = sorted((U*abs(i-j),V*abs(l-m)))
                if not d:
                    continue
                t,rem = divmod(d*e,N)
                assert rem == 0 and t >= 3
                k,rem = divmod(t*(d*d+e*e-t*t-1),N-1)
                assert rem == 0 and k>0
                assert (k-2*t*t)*(N-1) == t*((d-e)**2-(t-1)**2)
                assert t*(d*d-1)*(d*d-t*t) % (N-1) == 0
                assert (2*t*t*e-k*d)**2-(k*k-4*t**4)*d*d == 4*t**3*(t**3+t-k)
                total += 1
                resonance += k == 2*t*t
    print('Low-residue controls:',len(systems),'assignments,',total,'gap pairs,',resonance,'resonant')


if __name__ == '__main__':
    check_coefficients(45)
    check_binary_cutoff(1024)
    check_low_residues()
