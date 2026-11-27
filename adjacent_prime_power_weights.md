# Weights built from adjacent prime powers with an arbitrary multiplier

Let $p$ be prime, let $R=p^h-1>1$ be a prime power, and take integers $d,k\ge1$, $j\ge0$. At
$$
                     N=1+dR p^{hk+j}
$$
**Wasserman's three-row conjecture holds whenever**
$$
                     R>(k+3)d^4 p^{3j}2^k. \tag{1}
$$
There is no restriction on the number of prime divisors of $d$. With no lower bound on $R$, the result holds for **all $d\le512$, $1\le k\le6$ and $0\le j\le\min(h-1,2)$**.

For example take $R=2^{521}-1$. The theorem applies simultaneously to every $d\le2^{96}$, every $1\le k\le48$, and every $0\le j\le24$. Furthermore, each verified exponent $hk+j$ gives an entire arithmetic progression of good exponents at this same cofactor $dR$, by [transfer between separated exponents](digit_gap_transfer.md).

## A necessary rational root

More generally, suppose $R=q^\beta>1$, $p^h=cR+1$, for primes $p,q$ and a positive integer $c$. Put
$$
 W=p^j(cR+1)^k,\quad A=dp^j,\quad \alpha=A c^k,\quad L=A(c+1)^k.
$$
Write $\mu_p(R)=\min_{t\ge1}s_p(tR)$ for the minimum base-$p$ digit sum over positive multiples. Assume
$$
 s_p(dR)\le 2\mu_p(R),\qquad R>(k+3)d\alpha^2 L. \tag{2}
$$
If the conjecture fails at $N=1+dRW$, there is a positive integer $v$ satisfying
$$
                         A(v-c)^k=v^{k+1}. \tag{3}
$$

In a counterexample the [endpoint row](all_weights_endpoint_reduction.md) occurs. The other two rows must avoid $p$ and $q$ between them, since both primes divide $N-1$. They cannot assign both primes to one row. Such a row would have two positive parts divisible by $R$ and a third divisible by its assigned weight scale, a divisor of $N$ at least three and coprime to $p$. The first two parts have digit sum at least $\mu_p(R)$ each. The third cannot be a power of $p$, so its digit sum is at least two. Carry-free addition in base $p$ would now require $1+s_p(dR)\ge2+2\mu_p(R)$. These digit sums also apply when $p$ or $q$ divides $d$, because $p\mid W$ and the full assigned powers only increase.

Call the $p$-avoiding row's weight scale $U$ and the $q$-avoiding row's scale $V$, so $UV=N$. At least two parts of the first row are positive multiples of $WU$, and at least two of the second are positive multiples of $RV$. Each part in either row is at least its scale. Consequently
$$
                   U>2R,\qquad V>2W,\qquad 2U<dR.
$$
For the last inequality, $N\ge2WU+U$ and $U\ge3$.

Set $H(X)=1+AX(cX+1)^k$, so $H(R)=N$. Its coefficients above the constant are the nonnegative coefficients of $A(cX+1)^k$ and have sum $L<R$. Thus they are the radix-$R$ digits of $N$ above its units place. Addition free of base-$q$ carries is free of carries in radix $R=q^\beta$ as well. Choose a positive part divisible by $RV$ in the second row. It equals $RF(R)=R\ell V$, where $1\le\ell<d$, and the coefficients of $F$ are nonnegative and at most the respective coefficients of $A(cX+1)^k$. In particular $F(R)=\ell V>2W\ge2R^k$, so $\deg F=k$.

Let $v_0,w_0$ be the two highest coefficients of $F$. They satisfy $1\le v_0\le\alpha$, $0\le w_0\le k\alpha$. The two highest coefficients of $H$ are $\alpha$ and $b=kAc^{k-1}\le k\alpha$. Cancel the first two terms by forming
$$
 D(X)=\ell v_0^2 H(X)-\ell(v_0\alpha X+v_0 b-\alpha w_0)F(X).
$$
Then $\deg D\le k-1$ when $D\ne0$. Its sum of absolute coefficient values is bounded by
$$
 \ell v_0^2(1+L)+\ell(v_0\alpha+|v_0b-\alpha w_0|)L
 \le(k+3)d\alpha^2L<R.
$$
Yet $F(R)$ divides $D(R)$ since $F(R)=\ell V$ divides $\ell H(R)$. Now $|D(R)|<R^k<F(R)$, forcing $D(R)=0$. If $D$ were nonzero, its lowest nonzero coefficient would then be a multiple of $R$, exceeding the norm bound. Therefore $F$ divides $H$ over the rationals and their degrees differ by one. The rational root of $H$ is negative since $H$ has positive coefficients, and has numerator $-1$ in lowest terms since $H(0)=1$. Substituting $-1/v$ gives (3).

One can list the possible multipliers in (3) without factoring the weight. Put $g=\gcd(v,c)$ and $u=v/g$. Then $|u-c/g|^k$ divides $g$, by coprimality, and $A=g u^{k+1}/(u-c/g)^k$. In particular for squarefree $c$ and $k\ge2$, only $u=c/g+1$ or $u=c/g-1$ need be considered.

For $c=1$, (3) forces $v-1=1$ by coprimality and $A=2^{k+1}$. Thus $d$ would be a power of two. One of $R,R+1$ is even, so $N-1$ would have at most two distinct prime divisors, an [already proved case](endpoint_reduction.md). Also $\mu_p(R)=h(p-1)$ by folding blocks of $h$ digits. If $dR<p^{3h}-1$, subtract $dR$ from $p^{3h}-1$. The difference is a positive multiple of $R$, giving $s_p(dR)\le2h(p-1)$ by complementary digits. The bound (1) implies $d<R$ and hence this inequality. Thus (2) gives (1). [Digit-sum proof](sturdy_cofactors.md)

## Checking the small markers

For the finite parameter range above, $p^j\le4$. When $R$ is odd, $p=2$. When $R$ is even, $R$ is a power of two. If $h>1$ has an odd prime factor, factoring $p^h-1$ supplies a nontrivial odd geometric factor. Otherwise $h$ is even, and the factors $p^{h/2}-1,p^{h/2}+1$ are powers of two separated by two, hence 2 and 4. This gives $R=8,p=3,h=2$.

It suffices to list powers of two and their predecessors through $9\cdot512^4\cdot4^3\cdot2^6$, retaining pairs of prime powers. The list of $R$ is
$$
 2,\ 3,\ 4,\ 7,\ 8,\ 16,\ 31,\ 127,\ 256,\ 8191,\ 65536,\ 131071,\ 524287,\ 2147483647.
$$
For all but the last, use (1), the [one-billion verification](full_range_1b.md), or the two-factor predecessor criterion where applicable. If the digit-sum inequality fails, test the weight by the full endpoint CRT enumeration. Otherwise split $U,V$ as in the proof. The top digits of the first row give positive multiples $x,y$ of $U$ with $x,y,dR-x-y$ adding without a carry at $p$. The second row has normalized parts $Ra,Rb,U-R(a+b)$. Test the actual second row at $q$ and at every prime assigned from $V$. Any survivor goes through the full CRT enumeration. Of 82,944 parameter triples, 71,488 scale incidences pass the size test, 7,362 pass the first binary or odd-base digit split, and 331 distinct weights counted within their respective markers reach the full CRT check. None of its 1,219 size-feasible assignments has a row on the larger-modulus side.

For the last marker put $R=2^{31}-1$. Its digit-sum inequality holds since $dR<2^{93}-1$. The second radix-$R$ digit of $N$ is $A=d2^j<R$. Choose the smaller of the positive $a,b$ from its second row, naming it $a$. Both $Va$ and $Vb$ are nonzero modulo $R$, since $R$ is prime, $R\nmid V$, and $a,b<d<R$. Their least residues sum to at most $A$. Thus
$$
 1\le z=(Va\bmod R)<A,\qquad zU\equiv a\bmod R,\qquad 4a<d.
$$
Reduce $a,z$ by their gcd, which preserves the congruence. Writing $U=\ell R+(a/z\bmod R)$ gives $2a\le\ell<d/2$, also valid after reduction. The binary marker makes $U$ odd and its nonzero residue makes it coprime to $R$. For each of these scales and each $k$ a necessary condition is
$$
             A\equiv-\big(R(R+1)^k\big)^{-1}\pmod U.
$$
Here $A\le2048<U$, so the residue recovers $A$ exactly. The C++ verifier checks 10,197,384 scale candidates. Just 17 scale and weight incidences pass the modular and size checks. None permits two disjoint positive multiples of $U$ within the binary upper digits of $dR$.

The Python verifier uses python-flint, recomposes the factorizations, and checks primality of their bases. It verifies primality of the 521-bit example with FLINT and independently checks its Lucas–Lehmer residue. Run both verifiers from the repository.

~~~sh
c++ -std=c++17 -O3 check_marker_split_scales.cpp -o /tmp/wasserman-marker-scales
python3 check_adjacent_marker_weights.py --scale-checker /tmp/wasserman-marker-scales
~~~
