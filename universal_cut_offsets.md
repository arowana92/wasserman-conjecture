# The closest exceptional cuts to rational fractions of the weight

Let $M=N-1$, $C_a=\binom Na$ with $0<a<N$, and let $K_B$ be the coefficient of a positive trinomial row of weight $N$. A cut is universal at its weight if $\gcd(NM,C_a,K_B)>1$ for every such $B$. This property covers every trinomial obtained by refining the cut. By the [endpoint reduction](all_weights_endpoint_reduction.md), each cut of either nonendpoint row of a counterexample must fail universality.

**Across all weights**, the following are the sharp minima of the **nonzero** displacement $|ta-N|$ for a cut failing universality. For $t=2$ the reflected cut at the same weight is equivalent.

| $t$ | Minimum $|ta-N|$ | Unique weight $N$ | Cut $a$ | Number of unordered witnessing trinomials for this cut |
|---:|---:|---:|---:|---:|
| 2 | 4,616 | 5,752 | 568 | 5 |
| 3 | 4,048 | 5,752 | 568 | 5 |
| 4 | 113 | 21,809 | 5,424 | 1 |
| 5 | 2,912 | 5,752 | 568 | 5 |
| 6 | 2,344 | 5,752 | 568 | 5 |

In particular, at every odd weight all cuts with $|2a-N|<4616$ are universal. At every weight an arithmetic progression of parts with middle entry $b$ has $3b=N$; perturbing the two successive gaps so that their difference has magnitude from 1 through 4,047 guarantees a universal cut at $b$. Both statements allow an arbitrary opposite trinomial.

**At any fixed rational fraction $r/t$ in lowest terms**, an exceptional cut must satisfy $|ta-rN|\gg_{r,t} N^{2/7}$ outside the three exact offsets $0,-r,t-r$. For each fixed $L>0$, at most $O_{r,t,L,\varepsilon}(X^{2/7+\varepsilon})$ weights $N\le X$ can have an exception outside these exact offsets with $|ta-rN|\le L N^{2/7}$. There is also a growing common-divisor bound below.

**For exact fractions of the predecessor with reduced denominator at most seven**, there is a single failure up to the following limit. Write $N=tm+1$, $2\le t\le7$, $1\le m\le10^{1000}$, and $a=km$ with $0<k<t$ coprime to $t$. The only exception is $(t,k,m)=(7,5,8452)$. This includes all factorizations of $m$. At that weight the two witnessing trinomials, up to permutation, are $11833(1,1,3)$ and $11833(1,2,2)$. Reflection gives the equivalent cut $16905$ at $N=59165$.

## A divisor bound at one cut

Put $d=ta-rN$ with $0<r<t$. Let $D=\gcd(N,C_a,K_B)$ and $J=\gcd(M,C_a,K_B)$, so the common endpoint divisor is $G=DJ$. Use the [retained-power assignment](endpoint_pair_bound.md) for the binomial and the trinomial. Write $DUV=N$, $JST=M$, assigning $U,S$ to the cut. Then $U\mid a$, and at every retained full power in $S$ the index $a$ is zero or one. Consequently
$$ U\mid d,\qquad S\mid L_d=\operatorname{lcm}(|d+r|,|d+r-t|) $$
when $d+r$ and $d+r-t$ are nonzero. Indeed $N\equiv1$ at each such power, so it divides $d+r$ or $d+r-t$. The companion has parts $Vy_i$ with $\sum y_i=DU$. At a power in $T$ two $y_i$ vanish modulo that power, giving $27T^2\le (DU)^3$. Hence
$$27 M^2\le J^2D^3\gcd(|d|,N/D)^3\gcd(L_d,M/J)^2. \tag{1}$$
For $d\ne0$ this implies $G^3 |d|^3 L_d^2\ge27M^2$. At any fixed offset outside the three exact cases, the common divisor grows at least as a positive constant times $N^{2/3}$ for **every** companion. For $G=1$ the inequality gives the asserted $N^{2/7}$ scale.

For a fixed nonexact $d$, all possible weights admitting an exception can be searched finitely. Enumerate $U\mid |d|$, $U\ge3$, and $S\mid L_d$ coprime to $U$. Then
$$1\le T\le\sqrt{U^3/27},\qquad ST\equiv-1\pmod U,\qquad N=ST+1.$$
This leaves $O_{r,t,\varepsilon}(|d|^{1/2+\varepsilon})$ choices counting assignments, by the divisor bound and the progression step $U$. Thus across nonexact offsets of magnitude at most $H$, at most $O_{r,t,\varepsilon}(H^{3/2+\varepsilon})$ weights can occur, with no upper cutoff on $N$.

## Counting weights near the $2/7$ boundary

Take $G=1$ and $|d|\le L N^{2/7}$ outside the exact forms. It suffices to count large $N$. Since $M=ST$, $U\le |d|$, $T\ll U^{3/2}$, and $S\ll_{r,t} d^2$, the factors must have sizes
$$ U\asymp |d|\asymp N^{2/7},\qquad S\asymp U^2,\qquad T\asymp U^{3/2}. $$
In particular $h=d/U$ and $k=(d+r)(d+r-t)/S$ range over bounded sets of integers, with $k>0$ for large $N$. Put $R=r(t-r)>0$. From $ST\equiv-1\pmod U$ and $kS\equiv-R\pmod U$, the quantity $z=(RT-k)/U$ is a positive integer of order $\sqrt U$.

Group the powers in $T$ by their exceptional positions in $B$ with products $b_1,b_2,b_3$. Write $y_i=f_iT/b_i$. Each $f_i$ is a positive integer, and $\prod_i f_i\le U^3/(27T^2)$ bounds all three. Each $b_i>T/U$ since $y_i<U$, hence each is of order $\sqrt U$. The two identities for $T$ and $U$ give
$$ R b_1b_2b_3-z(f_1b_2b_3+f_2b_1b_3+f_3b_1b_2)=k.$$
Fix the bounded parameters and the two integers $z,b_1\ll X^{1/7}$. With $A=Rb_1-z f_1$, this equation gives $A>0$ and
$$ (A b_2-zb_1 f_2)(A b_3-zb_1 f_3)=kA+z^2b_1^2 f_2 f_3. $$
Both factors are positive, since, for example, multiplying the first by $b_3$ gives $k+zb_1 f_3 b_2$. The right side has $O_\varepsilon(X^\varepsilon)$ divisors. A factor pair determines $b_2,b_3$, then $T,U,d,S,N$. This proves the count.

## The exact offsets

An exceptional cut has $U\ge3$, since the positive normalized companion parts sum to $U$. If $r=1$, offset $-1$ is therefore impossible. At offset $t-1$ write $N=tm+1$, $a=m+1$; its other side is $(t-1)m$. Here $U\mid\gcd(m+1,t-1)$. If $t-1$ is 1 or 2 this already excludes the case. If $t-1=p$ is any prime, an exception would require $U=p$, and $p\mid m+1$. Put $j=v_p(m+1)$. The digit of $p m$ at place $j$ is $p-1$, because the low $j$ digits of $m$ are all $p-1$. The nonzero digit of $m+1$ there forces a carry in the very base assigned to the cut. The same argument works for $t-1=4$: then $U=4$ and $j=v_2(m+1)\ge2$, so $4m$ has a one at place $j$. This handles both exceptional predecessor offsets for the five minimum-displacement calculations. The case $ta=N$ is outside those minima.

For the full predecessor-fraction assertion, a cut $a=km$ at $N=tm+1$ requires $U\mid k$ since $\gcd(m,N)=1$. This eliminates $k\le2$, and the preceding paragraph and reflection cover $k=t-1$ for $t\le6$. The remaining reduced pairs are $(t,k)=(5,3),(7,3),(7,4),(7,5),(7,6)$.

At $(7,4)$ the scale must be four. Thus $m$ is odd and the two sides $4m,3m+1$ would have to add without a binary carry. This is impossible: for disjoint positive binary integers $x\ge y$, the number of ones of $x-y$ is at least that of $x$. Subtracting from $x$ each occupied bit of $x-y$ clears at most one current one, and all original ones must disappear to obtain $y$. Here $x-y=m-1$ has one fewer one than $x=4m$. This is the [binary difference observation](predecessor_five_seven.md#full-factor-seven-at-even-weights).

For the other four pairs set $p=k$ for $k=3,5$ and $p=3$ for $k=6$. Always $p\mid U$. For $k=3,5$ the companion is a multiple of a positive partition of 3 or 5 and must repeat a part, which makes its coefficient even. For $k=6$, either $U=3$ and the same argument applies, or $U=6$ assigns 2 to the cut. Thus the cut must avoid both 2 and $p$. It also must avoid $t$: a companion avoiding this predecessor prime would have two normalized parts divisible by $t$, but $U\le k<t$. At its full power in $M=tm$, the cut cannot be zero because $k$ is coprime to $t$. Its residue must be 1, forcing $km\equiv1\pmod t$ and $t\nmid m$. Additionally $tm+1\equiv0\pmod p$.

The resulting carry computations use the two summands $km,(t-k)m+1$. The table gives the complete results for $2\le m\le10^{1000}$. The necessary congruences also exclude $m=1$.

| $(t,k)$ | Required bases | $m$ modulo $tp$ | Interval nodes | Direct checks | Survivors |
|---|---|---|---:|---:|---|
| $(5,3)$ | $2,3,5$ | $7\pmod {15}$ | 118,280 | 31 | none |
| $(7,3)$ | $2,3,7$ | $5\pmod {21}$ | 148,952 | 39 | none |
| $(7,5)$ | $2,5,7$ | $17\pmod {35}$ | 561,596 | 54 | 8,452 |
| $(7,6)$ | $2,3,7$ | $20\pmod {21}$ | 496,488 | 91 | none |

The [digit-graph intersection algorithm](carry_intersections.md#exhaustive-search-by-intervals) tests the full three-base conditions on the indicated residue classes. These carry graphs arose previously for three-prime predecessors; the cut argument above makes them necessary for unrestricted $m$. At the one survivor $N=5\cdot11833$ and $U=5$, so the two positive partitions of five give every possible companion. Direct exact coefficients verify both gcds.

## Reproduction

The Python offset search deliberately uses divisors of the product in place of the lcm, a superset. It enumerates $T$ in its required progression, assigns every prime avoided by the cut to the cut itself, and tests the companion by CRT and Kummer digits. This produces respectively 5,294, 6,239, 9, 2,148 and 1,479 canonical companion tests at the five denominators. It also checks (1) against 367,414 exact coefficient comparisons at small weights, checks each boundary witness by full factorial coefficients, and runs the four shifted-fraction intersections.

The independent C++ program enumerates divisors of the lcm, generates $V$ via $UV\equiv1\pmod S$, and uses factorial floors. To check the companion it enumerates two positive normalized multiples of the largest full power in $T$, taking that multiple condition to be 1 when $T=1$. This covers all permutations and recovers every unordered witnessing row.

~~~sh
mkdir -p .work
c++ -O3 -std=c++17 check_universal_cut_offsets.cpp -o .work/universal_offsets
.work/universal_offsets
python3 check_universal_cut_offsets.py  # requires python-flint
~~~
