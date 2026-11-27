"""Check the two weight-factor obstructions for the residual quadratic shapes."""

from collections import Counter
from math import gcd
from check_seven_quadratic_suffixes import RESIDUAL
from check_endpoint_pair_bound import prime_powers


def census(a, b, m, n):
    # For the full power and cofactor f, the two numerators are 5 f^2 m,n.
    units = [f for f in range(a*b) if gcd(f,a*b) == 1]
    passed = [f for f in units if 2*(b*(5*f*f*m%a) + a*(5*f*f*n%b)) < a*b]
    return len(units), len(passed), 1 in passed


def exact_layers(bound):
    counts = Counter()
    for P in range(1,bound,2):
        T,W,H = 2*P+1,5*P+1,10*P+7
        for s,power in prime_powers(W):
            if s == 2:
                continue
            F = W//power
            E = power
            while E>F:
                k = F*F*(power//E)%5
                residue = F*P%E
                assert 5*residue == k*E-F
                assert (2*residue+F<E) == ((2*k-5)*E+3*F<0)
                counts['W'] += 1
                E //= s
        # The residue calculation is an integer identity also at composite R.
        for D,a,b,z,w in RESIDUAL:
            if H%D:
                continue
            R = H//D
            m,n = z//(D//a),w//(D//b)
            for s,power in prime_powers(T):
                if s==2:
                    continue
                B = T//power
                E = power
                while E>max(3,2*m,2*n)*B:
                    u = 5*B*B*(power//E)*m%a
                    v = 5*B*B*(power//E)*n%b
                    x,y,N = B*R*z, B*R*w, B*W
                    assert a*(x%E) == u*E+2*B*m and b*(y%E) == v*E+2*B*n
                    assert 2*(N%E) == E-3*B
                    full = b*u+a*v
                    lhs = 2*E*full + B*(4*b*m+4*a*n+3*a*b)
                    assert (x%E+y%E <= N%E) == (lhs <= a*b*E)
                    counts['T'] += 1
                    E //= s
    return counts


if __name__ == '__main__':
    expected = [(2760,368,False),(2332,108,False),(5984,304,True)]
    for (D,a,b,z,w),answer in zip(RESIDUAL,expected):
        m,n = z//(D//a), w//(D//b)
        assert gcd(z,w) in (2,6) and D%3 == 0
        actual = census(a,b,m,n)
        assert actual == answer
        print(f'D={D}: {actual[1]} / {actual[0]} cofactor classes pass; class 1: {actual[2]}.')
    counts = exact_layers(300000)
    print('Direct residue comparisons:', dict(counts))
    assert counts['T']>0 and counts['W']>0
