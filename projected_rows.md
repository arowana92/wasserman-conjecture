# Determinants for the two rows paired with the endpoint

Write $K_A=N!/(A_1!A_2!A_3!)$ for the coefficient of a positive row $A$ of sum $N$. Fix two such rows $A,C$, put $M=N-1$, and set
$$G=\gcd(NM,K_A,K_C).$$
The [endpoint reduction](all_weights_endpoint_reduction.md) says that proving $G>1$ for every pair at a weight proves the full conjecture there. The bounds below force two rows near prescribed rational proportions into exact integer forms. Applied to binary digits, they prove density one in every growing interval immediately above a fixed binary prefix.

## A bound for an arbitrary positive direction

Choose any positive integer vector $z$ of length three. Write $e_0=0$ and let $e_1,e_2,e_3$ be the three coordinate unit vectors. Define
$$L_A=\max_{0\le i<j\le3}|\det(z,A-e_i,A-e_j)|.$$

**Theorem.** Suppose none of the four vectors $A-e_i$ is proportional to $z$. Then
$$G^3 L_A^5\ \ge\ 27M^2. \tag{1}$$
If the three determinants involving $e_0$ are all nonzero, the stronger bound is $G^6 L_A^9\ge729M^4$. If $C$ also satisfies the nonproportionality hypothesis for a positive direction of its own, then
$$G\,L_A^2L_C^2\ \ge\ NM. \tag{2}$$

For example take $z=(1,1,1)$ and let $w_A=\max_i A_i-\min_i A_i$. The six absolute determinants are the three pairwise differences in absolute value and the three numbers $|M-3A_i|$. Thus $L_A\le2w_A+1$. Apart from rows $(t,t,t)$ and permutations of $(t+1,t,t)$, (1) gives a growing lower bound on $G$ whenever $w_A=o(N^{2/5})$. For three distinct parts the exponent improves to $4/9$. When both rows satisfy the hypothesis, (2) forces $(2w_A+1)(2w_C+1)\ge\sqrt{NM}$ in any counterexample containing the endpoint.

**Proof.** Put $D=\gcd(N,K_A,K_C)$ and $J=\gcd(M,K_A,K_C)$, so $G=DJ$. Use the [partial prime-power assignment](endpoint_pair_bound.md). For each full power $p^e$ in $N$ or $M$, remove its common portion $p^{\min(e,v_p(K_A),v_p(K_C))}$ and assign the remaining power to a row of smaller valuation. At that retained power the residues of the row sum to 0 for $N$ and to 1 for $M$. Indeed once the low residues exceed the target 0 or 1 at some level, they exceed it at all higher levels through $e$, each contributing at least one to the valuation.

Let the retained products assigned to $A$ be $U$ from $N$ and $S$ from $M$, and those assigned to $C$ be $V,T$, so $DUV=N$ and $JST=M$. Split $S=S_1S_2S_3$ by the position of the residue 1. With $R_0=U$ and $R_i=S_i$ for $1\le i\le3$, every entry of $A-e_i$ is divisible by $R_i$. In particular
$$R_i R_j\mid\det(z,A-e_i,A-e_j). \tag{3}$$
The parts of $C/V$ sum to $DU$, and each retained power of $T$ divides two of them. Their product is at least $T^2$. The arithmetic–geometric mean inequality gives
$$27M^2\le J^2 D^3 U^3 S^2. \tag{4}$$

Project along the line spanned by $z$ into a plane. No three images of the $e_i$ are collinear. Their four triangle determinants have absolute values $z_1,z_2,z_3,z_1+z_2+z_3$, all positive. The image of $A$ differs from each of these points by hypothesis. A zero determinant in (3) means that the line through two projected points contains the image of $A$. Two such pairs cannot share a vertex, since their lines meet only at that vertex. Hence one can pair all four indices using two nonzero determinants. Multiplication in (3) gives $US\le L_A^2$. Also $U\le L_A$ by any nonzero pair incident to 0. Thus $U^3 S^2=U(US)^2\le L_A^5$, proving (1) by (4).

If $U\le\sqrt{L_A}$ this last argument gives $U^3 S^2\le L_A^{9/2}$. Otherwise, when all three pairs incident to 0 are nonzero, $S_i\le L_A/U$ for each $i$, giving the same bound. For (2) apply the pairing argument to each row and multiply $US\le L_A^2$ and $VT\le L_C^2$.

## A gap around prescribed rational proportions

Fix positive primitive directions $z,z'$ for the two rows, let $B_z=\max_i z_i$, and use the respective determinants in $L_A,L_C$. If $G=1$, $M>B_z^3B_{z'}^3$, and
$$L_A^2<M/B_{z'}^3,\qquad L_C^2<M/B_z^3,$$
then, in some order, **one row is an exact integer multiple of its direction and the other is an integer multiple of its direction plus one coordinate unit vector**. In particular, if both rows have deviations $o(\sqrt N)$ from fixed positive rational proportions along a sequence of counterexamples, they must eventually have these exact forms. Indeed, deviations $o(\sqrt N)$ from $Nz/(\sum z_i)$ make the determinants $o(\sqrt N)$, by subtracting this multiple of $z$ and expanding in $e_i-e_j$.

To prove the gap, an unshifted row $A=t z$ has $\gcd(t,M)=1$. Each $S_i$ divides both unchanged direction coordinates at its zero positions, so $S\le B_z^3$. For a shifted row $A=t z+e_i$, $i>0$, we have $U\mid\gcd(z_j:j\ne i)$, since $U\mid A$ coordinatewise and $\gcd(t,U)=1$. Thus $US\le B_z^3 N$ in either exact form. A row outside these forms has $US\le L_A^2$. Apply the same bounds to $VT$ in the other row. Since $(US)(VT)=NM$, the displayed inequalities require both rows to have exact forms. Two unshifted rows would give $M=ST\le B_z^3B_{z'}^3$ and two shifted rows would give $N=UV\le B_zB_{z'}$. Both are excluded. For fixed directions the number of weights through $X$ satisfying these small-deviation hypotheses and admitting such a pair with $G=1$ is $O_{z,z'}(X^{\log_2\varphi})$, by the binary encoding below. At least one row must avoid 2, so two positive unshifted multiples of its scalar must have disjoint ones. The scalar determines the weight in each of the finitely many forms.

For the equal direction $(1,1,1)$ there is also a uniform statement at every $N\ge3$. If $G=1$, then $2\max(w_A,w_C)+1\ge(NM)^{1/4}$. This follows from (2) when neither row has an exact form. A shifted equal row has gcd one and cannot receive its required scale of at least 3. If one row is exactly equal, it receives no predecessor prime since its three residues agree. The other must receive $T=M$ and has $V\ge3$, so it cannot be equal either and its determinant bound gives $L_C^2\ge3M$ after choosing the names accordingly.

## Two blocks of digits in any prime base

Let $N=dQ+c$, where $Q=p^k$, $p$ is prime, $k,d\ge1$, $1\le c<Q$, and $p\mid N(N-1)$. Suppose
$$(dc)^5<27(N-1)^2. \tag{5}$$
In a counterexample the row $A$ avoiding $p$ can be taken different from the endpoint. Write $A_i=a_iQ+b_i$. By the carry criterion the $a_i$ split the digits of $d$ without a carry and the $b_i$ split those of $c$ without a carry. Every part of $A$ has fifth power greater than $27(N-1)^2$, by the [individual-part bound](endpoint_factors.md), so all three $a_i$ are positive. Write $a=gz$ with $g=\gcd_i a_i$, $z$ primitive, and $s=\sum z_i=d/g$. Then necessarily
$$ b=mz+e_i\quad\text{for some }0\le i\le3\text{ and integer }m\ge0. \tag{6}$$
In the shifted case $i>0$, the gcd of the two $z_j$ with $j\ne i$ must be at least three.

Indeed the determinants are unchanged when $A$ is replaced by $b$, since $A-b$ is a multiple of $z$. For the pairs incident to 0 their absolute values are $|z_j b_l-z_l b_j|$. For the other pairs they are $|z_l(c-1)-s b_l|$, where $l$ is the omitted nonzero index. All are at most $sc\le dc$. Equations (1) and (5) force $b-e_i$ to be proportional to $z$. Primitivity makes the scalar integral, and an unchanged nonnegative coordinate makes it nonnegative. In the shifted case write $A=t z+e_i$. The scale $U\ge3$ divides every entry, while $\gcd(U,t)=1$, proving the last restriction.

In particular $s$ divides $d$ and one of $c,c-1$. Also $p\nmid s$ because the lowest occupied digit of the upper parts sums without a carry. In base two the valuations of the three $z_i$ are distinct, so $s\ge1+2+4=7$. In base three at least one $z_i$ is divisible by 3, so $s\ge5$. In all other bases $s\ge3$. Thus (5) proves the full conjecture if the prime-to-$p$ parts of both $\gcd(d,c)$ and $\gcd(d,c-1)$ are smaller than the respective threshold. The earlier [binary block criterion](separated_binary_blocks.md) reaches exponent $2/9$ for $dc$ with its two gcd restrictions.

## Counting entire families of offsets

**Theorem.** Fix $d\ge1$. Among positive offsets $c\le H$, at most $O_d(H^{\log_2\varphi})$ can admit even one exceptional weight $N=d2^k+c$ with $c<2^k$ and (5), where $\varphi=(1+\sqrt5)/2$. Here $\log_2\varphi<0.695$. The bound counts an offset once over *all* such exponents. In particular, outside this set every sufficiently large member of the corresponding family satisfies the conjecture.

For each of the finitely many binary partitions of the upper block, (6) gives $c=sm$ or $sm+1$. Two unchanged lower parts $z_jm,z_lm$ must have disjoint binary ones. For any fixed positive integers $u,v$ the number of $m<2^n$ with $um$ and $vm$ disjoint is $O_{u,v}(\varphi^n)$. Here is an encoding proving this bound.

Cancel a common power of two and swap the multipliers so that they are $u'$ and $2^h v'$, with $u',v'$ odd. For $h=0$ the lowest one overlaps at every $m>0$. Otherwise read the $n$ bits of $m$ from low to high. Put $x=u'm$ and $y=v'm$. If bit $i-h$ of $y$ is 1, bit $i$ of $x$ is forced to 0 and we record 0. At every other step record bit $i$ of $y$, taking negative positions to be zero. This encoding is injective. Given the lower input bits, the next bit of either odd multiple determines the next input bit uniquely. Hence at a forced step use the known zero bit of $x$, and otherwise use the recorded bit of $y$. Every recorded 1 forces a recorded 0 exactly $h$ positions later if that position exists. Each of the $h$ residue classes of positions is therefore a string without adjacent ones. There are $F_{r+2}=O(\varphi^r)$ such strings of length $r$, by splitting on their last digit. Multiplying these bounds proves the count. For the pair $(1,2)$ the exponent is attained.

**An explicit prefix.** For $d=7$ the only positive upper parts, up to permutation, are $(1,2,4)$. Any two have gcd at most 2, so the shifted case in (6) is impossible. The remaining offsets are exactly $c=7m$ for which the ones of $m$ are at least three positions apart. Indeed $m,2m,4m$ must be disjoint. Conversely such an offset gives the binary carry-free proportional row $(2^k+m)(1,2,4)$ when $c<2^k$. It still has to meet the other prime conditions of a counterexample. The number of these offsets is $O(H^{\log_2\lambda})$ where $\lambda^3=\lambda^2+1$, $\lambda>1$, and $\log_2\lambda<0.552$. The recurrence for length-$n$ strings is $f(n)=f(n-1)+f(n-3)$ for $n\ge3$.

## Density in every growing interval after a fixed binary prefix

For each fixed $d\ge1$ and all integers $k\ge1$, $1\le H<2^k$, the number of exceptional weights in
$$d2^k<N\le d2^k+H$$
is $O_d(H^{15/16})$. In particular the exceptional proportion tends to zero even if $H$ grows arbitrarily slowly with $k$.

If $(dH)^5<27(d2^k)^2$, the offset theorem proves more. Otherwise $X\asymp_d 2^k$ satisfies $X\ll_d H^{5/2}$. We detail how to localize the [global exceptional-set estimate](endpoint_count.md). Split the weight interval at a dyadic boundary if needed so that $X/2\le N<X$. Enlarge the lower offsets to $0\le c<2^n$, where $n=\lceil\log_2(H+1)\rceil\le k$. Put $l=\lfloor10n/13\rfloor$. For the binary row in a counterexample its assignment has $K=US$. The small-modulus lemma there bounds the number with $K\le2^l$ by $X^{-1+o(1)}2^{5l/2}\ll_d H^{12/13+o(1)}$, using $X\ge H$ and $X\ll_d H^{5/2}$.

For $K>2^l$, fix the top $n-l$ bits of $c$. If they have $j$ ones and the bottom $l$ bits have $b$ ones, give each high one of $N$ weight 1 and each low one weight 2. There are $j+s_2(d)$ high ones, including the fixed prefix. Some part of the binary row has at most a third of the weighted total. With $A_0=3/2^{2/3}$, its possible values number at most $A_0^{j+s_2(d)}4^{b/3}(5/4)^l$, by the subset generating product evaluated at $1/2$. Each such part $w\ge3$ gives $X^{o(1)}$ possibilities for $U\mid w$ and $S\mid w(w-1)$. They fix $N$ modulo $K$ and hence its remaining low bits. We can also bound the number of weights by $\binom{l}{b}$. Take the geometric mean of these two estimates with exponents $3/5$ on the part count and $2/5$ on the binomial count. Use $\binom{l}{b}2^b\le3^l$ and sum over $b$ and all high strings. The resulting bound is
$$H^{o(1)}A_0^{3s_2(d)/5}(1+A_0^{3/5})^{n-l}(5/4)^{3l/5}3^{2l/5}.$$
Its exponent in $H$ is strictly below $15/16$ by the exact comparisons in the global estimate. In this application the fixed prefix contributes only a constant. This proves the claimed uniformity.

## Reproduction

Run `python3 check_projected_rows.py`. Exact binomial products at every weight from 3 through 36 check 89,671 ordered pairs of unordered rows. These give 266,666 nonproportional one-row bounds, including 7,495 which force a common factor by size, and 793,106 two-row bounds. The test includes 11,035 pairs with a partially retained prime power assigned to the first row. It verifies 52,388 binary encodings and their inverses, including every binary string with no adjacent ones for multipliers 1 and 2.

For $d=7$ a direct assignment of the nine low bits to parts agrees with the spacing test at every offset through 511. Among $1\le c\le10^6$ exactly 1,148 offsets remain after this test. Thus **for every $k\ge52$ the same 998,852 offsets give weights satisfying the conjecture**. The integer inequality $(7\cdot10^6)^5<27(7\cdot2^{52})^2$ works at the entire range of offsets.

The [classification of the shifted case](shifted_binary_prefixes.md) treats the possibility $b=mz+e_i$, $i>0$, in (6) and gives an exponent count for every fixed prefix.
