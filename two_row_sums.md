# A divisibility bound for two rows near a vertex

Let $A=(a,x,y)$ and $C=(b,z,w)$ be positive trinomial rows of weight $N\ge3$. Write $M=N-1$, $h=x+y$, $k=z+w$, $D=\gcd(N,K_A,K_C)$ and $J=\gcd(M,K_A,K_C)$, where $K_A,K_C$ are their coefficients. The following quantity is always a **positive integer**:
$$
\frac{D^3J^2\,xyzw\,(NM-ab)}{N^3M^2}. \tag{1}
$$
In particular, if $A,C$ occur with the endpoint in a counterexample, **every choice of two parts in each row** satisfies
$$
h^2k^2(h+k-1)>16N^2M^2,\qquad (hk)^{13}>2^{19}(NM)^{10}. \tag{2}
$$
Thus the product of the two sums grows at least as $N^{20/13}$, improving the [earlier $N^{3/2}$ bound](cross_prime_power_gaps.md). At least one of the two sums exceeds $8^{1/5}(NM)^{2/5}$. The [endpoint reduction](all_weights_endpoint_reduction.md) makes this a restriction on every possible counterexample, without assumptions on digit sums or factorization.

## The two integers giving the bound

Use the [retained powers](endpoint_pair_bound.md), with $DUV=N$ and $JST=M$. Each entry of $A$ is divisible by $U$, and each entry of $C$ by $V$. At each retained predecessor power in $S$ the residues of $A$ are $1,0,0$ in some order. The corresponding statement holds for $T$ and $C$.

First $N^2M\mid D^2Jxyzw$ since at least one of the two selected entries of the assigned row vanishes at each retained predecessor power. Second, $N^3M^2\mid D^3J^2 abxyzw$ since two of the three entries vanish there. For both divisibilities the weight factors and the predecessor factors are coprime. Multiply the first integer quotient by $DJ$ and subtract the second one. This gives (1), which is positive because $a,b\le N-2$.

We have $NM-ab=N(h+k-1)-hk$ and $xyzw\le h^2 k^2/16$. It follows that
$$D^3J^2h^2k^2(h+k-1)>16N^2M^2.$$
For the product bound interchange rows so that $h\ge k$. The last inequality gives $D^3J^2 h^3k^2>8N^2M^2$. At each retained power in $T$, either both $z/V,w/V$ vanish or exactly one vanishes and $k-1$ also vanishes. Thus
$$T^4\mid (zw/V^2)^3(k-1),\qquad S^2\mid axy/U^3.$$
Multiplying the first bound by the square of the second gives
$$
N^6 M^4\le D^6J^4(zw)^3(k-1)(axy)^2.
$$
Use $zw\le k^2/4$, $axy<Nh^2/4$ to get $D^6J^4 h^4k^7>2^{10}N^4M^4$. Multiply this inequality by the cube of the bound on $h^3k^2$. The result for arbitrary rows is
$$D^{15}J^{10}(hk)^{13}>2^{19}(NM)^{10}.$$
In a counterexample $D=J=1$, proving (2).

## Intervals at every degree from five upwards

Let $m\ge5$ be an integer and $P<Q$ powers of distinct primes, with $s=Q-P$ and $1\le s\le Q/(4m)$. The full conjecture holds at **every multiple of $P Q^{\lfloor m/2\rfloor}$** in the interval
$$
Q^m\le N\le Q^m+\frac{Q^{(m+3)/2}}{(ms)^{3/2}}. \tag{3}
$$

Write $h=N-P^m$, $k=N-Q^m$. We have $k\le Q^{m-1}$ and $P^m\ge(1-ms/Q)Q^m\ge3Q^m/4$, so both powers lie above $N/2$. A row avoiding the base prime of either power has a part at least that power. Its two other parts have sum at most the corresponding remainder.

We have $h\le (m+1)sQ^{m-1}$ and $h+k\le(m+2)sQ^{m-1}$. Consequently
$$h^2 k^2(h+k)\le \frac{(m+1)^2(m+2)}{m^3}Q^{4m}<3Q^{4m}<16N^2(N-1)^2.$$
By (2), in any endpoint counterexample the two base primes cannot be avoided by different rows. Both divide $N$, so both must be avoided by the same row.

Put $j=\lfloor m/2\rfloor$. Every entry of this row is divisible by $PQ^j$. The two smaller entries are $PQ^j x,PQ^j y$ with $x,y>0$, $x+y\le k/(PQ^j)$. In base equal to the prime below $Q$, division by $Q^j$ preserves carry freedom. But $P\equiv-s\bmod Q$, so the residues of $Px,Py$ modulo $Q$ have sum at least $2Q-s(x+y)$. This is at least $Q$, because $(m+3)/2\le j+2$, $P\ge Q/2$, and (3) gives $sk\le P Q^{j+1}$. The resulting carry is a contradiction.

There are infinitely many weights in (3) even if both leading remainders are required to be large. Fix $m$, take consecutive primes $P,Q$ in $[Y,2Y]$, and restrict to $s\le C\log Y$. By the prime number theorem there are $\gg Y/\log Y$ such pairs for a sufficiently large absolute $C$, since the sum of all consecutive gaps in the interval is at most $Y$. In the upper half of (3), each pair gives
$$\gg_m Y^{(m+3)/2-\lfloor m/2\rfloor-1}/(\log Y)^{3/2}$$
multiples. The exponent is one at odd $m$ and one half at even $m$, so the interval eventually contains this many integers in the required progression. Each weight is $O_m(Y^m)$ and hence has $O_m(1)$ prime factors in $[Y,2Y]$, bounding repetitions. In particular at degree five there are $\gg X^{2/5}/(\log X)^{5/2}$ such weights through $X$. At those weights **both leading remainders are at least a constant times $N^{4/5}/(\log N)^{3/2}$**. The earlier split-row product test alone does not cover this scale.

For a numerical example take $P=1{,}000{,}000{,}007$, $Q=1{,}000{,}000{,}009$ and degree five. The largest allowed multiple in (3) is
$$N=1000000045031623586790576688545947310930172393.$$
The verifier checks both primes, both actual remainders, the modulus-$Q$ obstruction, and the strict inequality $hk>4N\sqrt{N-1}$ at this weight.

Run `python3 check_two_row_sums.py`. Using exact binomial products, the program checks (1) and the product bound on 1,228,068 choices from 136,452 row pairs through weight 45. In 2,756 choices the old four-part size bound permits gcd one and the size bound from (1) excludes it. The separate interval check tests 73,995 upper endpoints across degrees 5 through 17.
