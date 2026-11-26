# Seven digit units on three positive powers

**Theorem.** For every prime $p$ and all six positive exponents $e_1,\ldots,e_6$ taking at most three distinct values, Wasserman's three-row conjecture holds at
$$N=1+\sum_{i=1}^6p^{e_i}.$$
For example, $1+2p^a+p^b+3p^c$ works in every ordering of $a,b,c$, as do $1+4p^a+p^b+p^c$ and $1+2(p^a+p^b+p^c)$. There is no upper bound on the base, exponents or gaps.

If the formal units carry and $p\ge3$, at most five units remain, so the [five-unit theorem](six_digit_gaps.md) applies. At $p=2$, every partition of six into at most three multiplicities has sum of binary digit sums at most four; hence $N$ has at most five ones. In what follows $p$ is odd and the seven units add freely.

## Two pairs at an odd predecessor prime

More generally let $p$ be odd, $p\mid N-1$ and $s_p(N)=7$, without any restriction on the number of occupied positions. In a counterexample the [endpoint theorem](all_weights_endpoint_reduction.md) supplies a row with coefficient $N(N-1)$, so a nonendpoint row $A$ must avoid $p$. In the [factor assignment](endpoint_reduction.md) assign $p$ to $A$ and let $U,V$ be the weight factors for $A,C$ respectively, where $C$ is the last row. Then $UV=N$, $U,V\ge3$, $U$ divides every part of $A$, and $p\nmid U$. No part of $A$ can consist of a single digit unit. Kummer's criterion therefore partitions the seven units into two pairs and one triple.

Write the paired parts as $p^i(1+p^d)$ and $p^j(1+p^e)$. Both $d,e$ are positive since a repeated power would force $U\mid2$. For $g=\gcd(d,e)$ both $d/g,e/g$ are odd and $p^g\equiv-1\bmod U$. Indeed the unit $p^g$ has square one modulo $U$ by Bézout applied to $2d/g,2e/g$. An even quotient would then give $1=-1\bmod U$, impossible. The paired parts have the same 2-adic valuation, since $(x^t+1)/(x+1)$ is odd at odd $x,t$. Their addition has a binary carry. Consequently $C$ must avoid 2. Since $N$ is odd, $V$ is odd and the three parts of $C/V$ have distinct 2-adic valuations. Their sum is $U$, whence **$U\ge7$**. In fact $U\ge1+3\cdot2^{v_2(N-1)}$, because the two even parts of $C$ have distinct valuations at least $v_2(N-1)$.

There is an exponent-only obstruction for any such seven-unit weight. Let $i_1\le i_2\le i_3$ be the triple's exponents and $t=\gcd(g,i_2-i_1,i_3-i_1)$. Then **$g/t\ge3$**. Put $y=p^t$. If $g/t=1$, reducing the triple after removing $p^{i_1}$ gives a sum of three signs modulo $U$, so $U\le3$. If $g/t=2$, use $y^2=-1\bmod U$ to obtain $a+b y=0\bmod U$ with $|a|+|b|\le3$ odd. If one coefficient vanishes the other has absolute value at most three. Otherwise $U\mid a^2+b^2\le5$. Both alternatives contradict $U\ge7$.

## The inner units

Write $W=p^{\min e_i}$ and $m=(N-1)/W$. The row $A$ has parts $(Wx,Wy,Wz+1)$ with $x+y+z=m$. All three are positive since $U$ cannot divide 1. Neither $x$ nor $y$ can use a single position: its multiplicity is at most six, while $U\ge7$ and $U$ is coprime to $p$. Split $m=ST$ by assigning the full prime powers in $S$ to $A$ and those in $T$ to $C$. Here $m$ is coprime to $p$, and $S$ is odd since $A$ carries at 2. Each power in $S$ divides $x$ or $y$, because two actual parts vanish there. Two parts of $C/V$ are divisible by every power of $T$. Their product, with the third normalized part, is at most $(U/3)^3$. Thus
$$27m^2\le U^3 S^2,\qquad U\mid\gcd(x,y),\qquad S\mid\operatorname{lcm}(\gcd(m,x),\gcd(m,y)). \tag{1}$$

One occupied inner position is already impossible. With two positions $1,J=p^n$, each of $x,y$ must take both, leaving at least one unit for $z$. The only possible unequal pair, up to interchanging the levels, is $J+1,2J+1$, whose gcd is 1 (the other ordering gives $J+1,J+2$). For $x=y=J+1$ we have $z=J+1,2$ or $2J$. The first forces $U\mid1$ from $Wz+1$. In the other two cases $S=1$ since $\gcd(m,J+1)\mid2$. Also $U\mid p^n+1$ and $U\mid2p^s+1$ where $s$ is the exponent of $W$ or $WJ$. Raising the latter congruence $2p^s=-1$ to the $n$-th power gives $U\le2^n+1$. Then $U^3<27p^{2n}\le27m^2$ for $p\ge3$.

For three occupied inner positions denote their values by $L,K,A_0$ in the order of multiplicity specified below. They are a permutation of $1,B^\beta,B^\alpha$, where $B$ is a power of $p$ and $\alpha>\beta\ge1$ are coprime. With multiplicities $3,1,2$ the nine possible $(x,y,z)$ are
$$
\begin{array}{c|ccc}
 &x&y&z\\\hline
1&K+A_0&L+A_0&2L\\
2&K+A_0&2L+A_0&L\\
3&L+A_0&L+A_0&L+K\\
4&L+A_0&L+K&L+A_0\\
5&L+A_0&L+K+A_0&L\\
6&L+A_0&2L+A_0&K\\
7&L+A_0&2L+K&A_0\\
8&L+2A_0&L+K&L\\
9&L+K&2L+A_0&A_0
\end{array}
$$
Cases 4, 5 and 6 respectively give $U\mid1$ from $x=z$, $U\mid K$ from $y-x$, and $U\mid A_0$ from $2x-y$. Case 3 is handled separately below. With multiplicities $4,1,1$ the two unshifted parts must take different singletons, so up to symmetry $(x,y)=(L+K,L+A_0)$ or $(L+K,2L+A_0)$.

With multiplicities $2,2,2$, choose a two-unit unshifted part, say $x=L+K$. A full power of $S$ cannot divide $x$, since $m=2(x+A_0)$ and $\gcd(S,2p)=1$. There are four cases for $y$ up to symmetry. If $y=L+K$ then $S=1$, $z=2A_0$ and $U\mid p^n+1$ for the gap $n$ between $L,K$. As in the two-position proof $U\le2^n+1$ from $U\mid2WA_0+1$. If $y=L+A_0$, again $S=1$ and $U\le B+1$, by Bézout applied to the two paired gaps, whose gcd is the exponent of $B$ in base $p$. Since $m\ge B^\alpha$ and $\alpha\ge2$, (1) fails. If $y=L+K+A_0$, then $U\mid A_0$. This leaves $y=K+2A_0$.

## Eight allocations bounded by two power congruences

The following table bounds $U$ and the parts of $S$ dividing $x$ or $y$. At each assigned modulus the coordinates $(L,K,A_0)$ must be proportional to the indicated vector, by a unit scalar. The symbol $-$ denotes modulus 1. The predecessor powers can be allocated to either eligible column when both apply.
$$
\begin{array}{c|cc|ccc}
(\#L,\#K,\#A_0)&x&y&U&S_x&S_y\\\hline
(3,1,2)&K+A_0&L+A_0&(-1,-1,1)&(-1,-3,3)&(-1,1,1)\\
 &K+A_0&2L+A_0&(-1,-2,2)&(-1,-3,3)&(-1,-1,2)\\
 &L+A_0&2L+K&(-1,2,1)&(-1,1,1)&(-2,4,1)\\
 &L+2A_0&L+K&(-2,2,1)&(-2,4,1)&(-1,1,1)\\
 &L+K&2L+A_0&(-1,1,2)&(-1,1,1)&(-1,-1,2)\\
(4,1,1)&L+K&L+A_0&(1,-1,-1)&(1,-1,-3)&(1,-3,-1)\\
 &L+K&2L+A_0&(1,-1,-2)&(1,-1,-3)&(1,-2,-2)\\
(2,2,2)&L+K&K+2A_0&(2,-2,1)&-&(1,-2,1)
\end{array}
$$
For the first data column set $x=y=0$ modulo the modulus, and for each predecessor column set $m=0$ and the corresponding part to zero. These linear equations give the vectors by substitution. All three levels are units at these moduli. Division by 2 is permitted because the moduli are odd. A factor of 3 in the vector is also a unit: the corresponding equation would otherwise set a level to zero. The first two multiplicity patterns require $p\ge5$ in the carry-free case; the last uses only coefficients prime to every odd $p$.

For a vector rearranged as $(c_0,c_1,c_2)$ in ascending order of the exponents, its modulus $Z$ divides
$$c_2^i c_0^{j-i} B^{(j\beta-i\alpha)_+}
       -c_1^j B^{(i\alpha-j\beta)_+}\qquad(0\le i\le j,\ j\ge1). \tag{2}$$
Indeed $B^\alpha=c_2/c_0$, $B^\beta=c_1/c_0$ modulo $Z$. The two coefficient magnitudes are at most $4^j$.

First take $(i,j)=(\beta,\alpha)$. Each nonzero difference bounds its modulus by $K_0=2\cdot4^\alpha$. Equality of the two terms implies $(c_2/c_0)^\beta=(c_1/c_0)^\alpha$. By Bézout there is a rational $\lambda$ with $c_2/c_0=\lambda^\alpha$ and $c_1/c_0=\lambda^\beta$, and $B=\lambda$ modulo $Z$. If $|\lambda|\ne1$, the three coordinate magnitudes must be distinct. The only such vector is $(-2,4,1)$, which requires the exponent at $L$ to be the midpoint between the exponents at $A_0,K$. Thus $KA_0=L^2$. In its two rows let $u=L/A_0\bmod U$. Either $u+1=2u+u^2=0$, or $u+2=u+u^2=0$. These force $U\le2$, impossible.

Otherwise $\lambda=-1$ since the signs are mixed. The modulus divides $B+1$. At most one of the three columns in a row can give such a zero: only the first row has two vectors with all magnitudes equal, and they prescribe different relative signs. The signs of the powers of $-1$ at the three fixed positions permit only one prescription up to a global sign. Accordingly $U^3S^2$ is bounded by one of
$$K_0^7,\qquad (2B)^3K_0^4,\qquad (2B)^2K_0^5.$$
At $B\ge256$ and $\alpha\ge4$, division by $B^{2\alpha}$ bounds these expressions by $2^{7-2\alpha}$, $2^{31-8\alpha}$ and $2^{23-6\alpha}$ respectively. Each is less than 27. For $(\alpha,\beta)=(2,1),(3,1),(3,2)$ the verifier evaluates the exact coefficient differences for all eight rows and all six orders. There are 140 inequalities and four exclusions from $KA_0=L^2$. With a zero replaced by $B+1$, the resulting ratio to $B^{2\alpha}$ decreases for $B\ge256$, because the exponent on $B+1$ is at most three. All 140 pass at 256.

For $\alpha>5$, two of the six fractional parts of $r\beta/\alpha$, $0\le r\le5$, lie in the same fifth of $[0,1)$. Their difference gives $0\le i\le j\le5$, $j\ge1$, with $|i\alpha-j\beta|\le\alpha/5$. This difference of exponents is nonzero by coprimality. In (2) the two terms have different $p$-adic valuations, so each modulus is at most $2^{11}(B^\alpha)^{1/5}$. Thus $U^3S^2\le2^{77}(B^\alpha)^{7/5}<27B^{2\alpha}$ when $B^\alpha\ge2^{128}$.

It remains to take $B<256$, $B^\alpha<2^{128}$, including the remaining $\alpha\le5$. We check even the larger set of odd integers $B\ge3$ in this range, allowing multiples of three in the last pattern and requiring $\gcd(B,6)=1$ in the other two. For each assignment remove all factors dividing $Bm$ from $\gcd(x,y)$, giving an upper bound $U_0$. Take the odd part of the lcm in (1) as $S_0$. Of 600,330 assignments, 58,001 have $U_0\ge7$; each fails $27m^2\le U_0^3 S_0^2$. An independent calculation takes the odd gcd of the three two-by-two minors of the vector of actual levels and each table vector. Using these three upper bounds gives the same counts and exclusion.

## The repeated part at multiplicities three, one and two

In case 3 the two unshifted parts are both $W(L+A_0)$. Each is a two-unit pair, while the triple contains 1. Their two gaps are equal, which contradicts the [general seven-unit gap inequality](seven_pair_gaps.md). This completes the last allocation.

## Reproduction

With python-flint installed, run `python3 check_seven_digit_units.py`. The program obtains all unshifted allocations by assigning the six labelled inner units to actual positions. It runs the two finite exclusions and independently computes the polynomial gcd for 30,856 exponent triples, including 378 nonconstant gcds. Numerical checks at three primes also verify (3) in the distinct-residue cases. At 728 sampled weight and base pairs a direct enumeration of all 33,326 positive carry-free rows verifies the scale or size obstruction, including weights with coincident levels.
