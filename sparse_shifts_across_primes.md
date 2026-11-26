# Counting large shifts across all prime bases

Fix $s\ge7$. Consider a tuple of $s-1$ positive integer positions $e_1\le\cdots\le e_{s-1}\le L$. Count a tuple once if, for **some** prime $p$ and positive integer $b$, the weight
$$N=1+\sum_{i=1}^{s-1}p^{b e_i}$$
has base-$p$ digit sum $s$ and could fail Wasserman's conjecture. Put $t=s-4$. The number of such tuples is at most
$$ O_s\big((1+L)^{t+1-1/(t+1)}\big). \tag{1}$$
The implied constant is effective. There is no restriction on the ratio of the smallest position to the span. For instance, eight digit units give $O(L^{24/5})$ tuples, out of order $L^7$ possible profiles. The [fixed-base count](effective_sparse_positions.md) is stronger for each individual prime. Here the same exceptional count covers the union over primes and overall dilations.

**Seven units.** The union in this case has the stronger bound $O((1+L)^3)$, even if the shift and inner dilation are independent. In detail consider $N=1+p^a\sum_{i=1}^6 p^{b h_i}$ for fixed $0=h_1\le\cdots\le h_6=H$. In any interval of length $H/9$, at most **4,815** distinct positive shifts $a$ can fail at any prime $p$ or dilation $b$ with digit sum seven. For instance, if $H\ge43{,}335$, every 4,816 consecutive shifts include one which works simultaneously for all such $p,b$.

## Reducing to one equation in a rectangular box

Write $P=p^b$, $a=e_1$, $h_i=e_i-a$ and $H=h_{s-1}$. The effective common-base gcd lemma and the restriction argument in the [sparse-position theorem](effective_sparse_positions.md#an-elementary-common-base-gcd-bound) apply with $P$ in place of $p$: all measured divisors are prime to $p$, $\mathbb Z[1/P]=\mathbb Z[1/p]$, and the proof uses integral powers of a base at least two. Recall two consequences of the theorem and its proof.

1. The $s-2$ inner coordinates lie in a fixed effective finite union of affine spaces of dimension at most $t$, even when inner dilation is allowed.
2. A counterexample with $H>0$ forces a divisor exceeding $P^{H/10}$ among the gcds of the three disjoint block equations in the [full-position argument](effective_sparse_positions.md#combining-primes-when-the-shift-is-at-most-a-fixed-multiple-of-the-span). Each eligible system has no common nonunit factor on a connected torus translate of codimension at most two containing a positive point. In the three-pair case the possibly degenerate system with right side $(0,0,0)$ is omitted. Its assigned factor is already at most $2m^{1/3}$, where $m=\sum P^{h_i}$, leaving one of the other systems with a large divisor.

An inner space of dimension less than $t$ contributes $O_s((1+L)^t)$ tuples after including the shift. On each remaining space, take an integral origin $h_0$ and a basis of its saturated lattice, parametrizing its integer points by $h=h_0+M w$, $w\in\mathbb Z^t$. Empty integer intersections may be discarded. In a range $R\le H<2R$, $R\ge1$, we have $\|w\|_\infty\ll_s R$ by a fixed left inverse. The map from $(a,w)$ to the full exponents has saturated image lattice: the coordinates $(a,h)$ come from an invertible integer change, and the shift is free. Its torus image is a positive translate of codimension two. Thus the restricted polynomials are collectively coprime. Clearing powers of $P$ and monomials, and choosing two coprime linear combinations as in the restriction lemma, gives bounded degrees and coefficient norms at most $P^{B_s}$ with effective constants independent of $p,b,a,H$.

For $a\le 2H$ the full-position theorem already gives $O_s((1+L)^t)$, including overall dilations by the same common-base argument. For $a>2H$ use dyadic intervals $A\le a<2A$, $R\le H<2R$, so $A\ge R$ and $A\le K=\|(a,w)\|_\infty\ll_s A$. Take $\epsilon=c_s R/A$ with a small effective positive rational $c_s$ so that $\epsilon K\le H/10$. Apply the explicit hyperplane lemma in $t+1$ variables to the large divisor. Its approximation denominator is at most $\lceil C_s A/R\rceil^t$. Hence either
$$ A\ll_s (A/R)^{t+1} \tag{2}$$
or there exist integers $u$ and $v\in\mathbb Z^t$, not both zero, such that
$$ |u a+v\cdot w|\le M_0,\qquad |u|,\ \|v\|_\infty\le M_0,
\qquad M_0=\lceil C'_s(A/R)^t\rceil. \tag{3}$$
Indeed the small-height alternative in the lemma has one extra factor $1/\epsilon$ compared with the bounds for the normal and offset. All the coefficient bounds used for (2) and (3) are uniform on the finite collection of spaces and block allocations.

## Counting the integer points

Let $R_*=A^{1-1/(t(t+1))}$. For $R\le R_*$ sum the trivial box bound $O_s(A R^t)$ over dyadic $R$. This gives $O_s(A^{t+1-1/(t+1)})$.

For $R>R_*$ one has $(A/R)^{t+1}<A^{1/t}$, so (2) is impossible once $A$ exceeds an effective constant. Also $M_0\ll_s A^{1/(t+1)}$, so $M_0\le R$ and $M_0<A$ for large $A$. We count (3) directly. Write $j=\|v\|_\infty$. When $v=0$ no nonzero $u$ works, as $a\ge A>M_0$. For $1\le j\le M_0$ there are $O_t(j^{t-1})$ vectors $v$ of norm $j$. Fix a coordinate attaining that norm. Once $a$ and the other $t-1$ coordinates of $w$ are fixed, there are $O(M_0/j)$ choices for the remaining coordinate in (3).

If $u=0$, this bounds the count by $O_s(A R^{t-1}M_0^t)$ after summing on $j$, since $t\ge3$. For nonzero $u$ a necessary condition is $|u|A\le O_s(jR)+M_0$. Since $M_0\le R$, there are $O_s(jR/A)$ nonzero integer choices for $u$, also valid when there are none. Summing these choices gives $O_s(R^t M_0^{t+1})$. This separates the cost of the long coordinate from the cost of the $t$ short coordinates.

On substituting $M_0\ll_s(A/R)^t$, both bounds decrease geometrically along dyadic $R>R_*$. At $R_*$ the second bound is $O_s(A R_*^t)$. The first is also at most this since $(A/R_*)^{t^2}=A^{t/(t+1)}\le R_*$. Summing in $R$, then in dyadic $A\le L$, proves (1). Bounded $A$ and $H=0$ add $O_s(1+L)$.

## The seven-unit union

Only the 26 effective bases $p^b<64$ can contribute, by the [large-base theorem](seven_refined_cutoffs.md). A genuine seven-unit binary expansion uses six distinct positive positions. For each remaining $p,b$, the original positions $a+b h_i$ are at most $O(L)$ if $a,H\le L$. Since the normalized list starts at zero, the map from $(a,h)$ to these positions is injective. Summing the fixed-prime $O(L^3)$ bound over these finitely many choices gives the stated combined count.

For fixed inner positions and an allocation, the [four-shift inequality](sparse_prime_counts.md#an-explicit-shift-bound-for-seven-units) permits at most three shifts in any interval of length $bH/9$, hence also of length $H/9$. At $(p,b)=(2,1)$ there are at most 105 allocations. For every other remaining base the unit at 1 belongs to a pair, by the [odd-base theorem](seven_large_bases.md) or [binary alignment](seven_binary_alignment.md). There are at most $6\binom52=60$ corresponding labelled allocations. The union bound is $3(105+25\cdot60)=4{,}815$. If $H=0$ there is at most one integer in an interval of length zero, so the spacing assertion also holds there.
