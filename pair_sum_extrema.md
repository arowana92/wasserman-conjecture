# Bounds for the sum of two parts near a prime power

Write $M=N-1$. Every counterexample includes the row $(1,1,N-2)$ by the [endpoint reduction](all_weights_endpoint_reduction.md). In either other row let $h$ be the sum of **any** two parts. The following are necessary conditions.

* If $M$ has at most four distinct prime factors, then $h^2(h-1)>16M^2$. With at most three distinct prime factors the stronger bound is $h^2(h-1)\ge36M^2$. At even $N$ the constants increase to 32 and 128 respectively.
* At every even $N$, $h^6(h-1)\ge32\cdot6^6M^4$. When $h$ is odd the factor 32 can be raised to 512. At $8\mid N$ it can be raised to 96 if the selected row has odd coefficient. The previous factor at all weights was one.
* At arbitrary $N$, set $H_M=6^{6/7}M^{4/7}$. Along any unbounded sequence of counterexamples,
  $$\liminf\frac{h-H_M}{H_M^{3/4}}\ \ge\ \frac{3\sqrt6}{7}.$$
  At the scale $h\le C N^{4/7}$ there are at most $O_{C,\epsilon}(X^{2/7+\epsilon})$ possible exceptional weights through $X$, for each fixed $C>0$ and $\epsilon>0$.

If $Q=p^k>N/2$ and $p\mid NM$, a row avoiding $p$ must place $Q$ in one part, and its two other parts sum to at most $N-Q$. It is one of the two nonendpoint rows. Thus each lower bound on $h$ gives a sufficient condition for the full conjecture in terms of the offset $N-Q$. In particular, for $\omega(N-1)\le4$ it suffices that
$$(N-Q)^2(N-Q-1)\le16(N-1)^2.$$
This includes weights with many digits in base $p$.

For an example take
$$N=1427247692705961760186075461832537708033263872,\qquad Q=2^{150}.$$
Here $M=1000000007\cdot2147483647\cdot4294967291\cdot154742504079676789$. The offset has 50 binary ones and exceeds $Q^{2/3}$. It meets the four-factor test. The verifier checks the four primes recursively and verifies the inequalities with integers.

For each sufficiently large power of two $Q$, the interval $Q+Q^{2/3}\le N\le Q+2Q^{2/3}$ contains $\gg Q^{2/3}/\log Q$ such weights with **exactly four** predecessor primes. Fix any three distinct odd primes and call their product $m$. Set $N=1+mr$ for primes $r$ in the corresponding interval. The short-interval prime number theorem supplies the count since its length is a fixed multiple of $(Q/m)^{2/3}$; the fourth prime exceeds the fixed three. The displayed test holds throughout the interval. [Guth and Maynard, Corollary 1.3](https://arxiv.org/abs/2405.20552)

## Factoring the two-part inequality

For any two rows $A,C$ of weight $N$, put $D=\gcd(N,K_A,K_C)$ and $J=\gcd(M,K_A,K_C)$. Assign the [retained prime powers](endpoint_pair_bound.md#retaining-partial-prime-powers) as $DUV=N$, $JST=M$, where the parts of $A$ are multiples of $U$, and those of $C$ are multiples of $V$. At a power assigned from $M$, the assigned row has one residue 1 and two residues 0. Write $S=dew$, placing each assigned power in the bin for its exceptional position, in order. The first two actual parts, denoted $x,y$, and their sum $h$ satisfy
$$x=Uw a e,\qquad y=Uw b d,\qquad h-1=cde$$
with $a,b,c$ positive integers. Split $T=t_1t_2t_3$ by the companion's exceptional positions, and write $C_i/V=f_iT/t_i$, $\ell=f_1 f_2 f_3$. Its three normalized entries sum to $DU$, so $27\ell T^2\le(DU)^3$. Put $P=w^2(ab)^3 c\ell^2$. Direct cancellation gives
$$
\frac{6^6 M^4 P}{D^6J^4 h^6(h-1)}
=\left(\frac{4xy}{h^2}\right)^3
 \left(\frac{27\ell T^2}{(DU)^3}\right)^2\le1. \tag{1}
$$
Thus (1) bounds a common endpoint divisor even for two arbitrary rows.

If $DJ$ is odd, assign 2 to a row which avoids it. At even $N$ all predecessor bins are odd. If $A$ takes 2, $U$ is even and $DU=\sum_i f_i T/t_i$ is even, so $\ell$ is even. The first two parts have different 2-adic valuations by carry freedom; therefore $ab$ is even and $P\ge32$. If $C$ takes 2, its three different valuations make $8\mid\ell$ and $P\ge64$. When $h$ is odd only the latter case can occur, and $ab$ must be even as well, giving 512. For odd $N$ and even $h$, $P\ge4$: if 2 is in $S$ it lies in $w$ since $de\mid h-1$. If it is in $T$, the two companion positions outside its bin have different 2-adic valuations but the same contribution from $T/t_i$, so $\ell$ is even.

For the remaining assertions use $D=J=1$, which holds for the two nonendpoint rows in a counterexample. We always have $U,V\ge3$.

## Above the two leading digits

Suppose $N=L+B+z$ where $L>B$ are powers of a single prime $p\mid N(N-1)$, $0\le z<B$, and $B+z<L$. In a counterexample the row avoiding $p$ places $L$ in one part. Write $x,y$ for the other two. Then $x+y\le B+z$ and $xy\le Bz$. Indeed at most one of $x,y$ receives the digit unit $B$. If one receives it, that part is at least half the sum bound, and the product is maximized at $(B,z)$. If neither receives it, $xy\le z^2/4$. Thus the original two-part inequality proves the full conjecture whenever
$$ (Bz)^3(B+z-1)<729(N-1)^4. $$
For example, **every integer** $2^{5n}+2^{3n}+z$ with $n\ge1$ and $0\le z\le2^{8n/3}$ satisfies the conjecture. The interval starts at an offset of order $N^{3/5}$ above the leading power and has length of order $N^{8/15}$. The estimate uses $B^4 z^3\le L^4$.

This also enlarges the initial interval above **every** power of two by a multiplicative factor. Write $k=7m+i$, $m\ge1$, $0\le i<7$. Every weight
$$2^k\le N\le2^k+(A_i/1000)\,2^{4m}$$
satisfies the conjecture, where in order of $i$ the $A_i$ are
$$5291,\quad6968,\quad11193,\quad15257,\quad23875,\quad40331,\quad51368.$$
In fact put $h_i=(6^6 2^{4i})^{1/7}$, $b_i=2^{\lfloor\log_2 h_i\rfloor}$. One may replace $A_i/1000$ by the root $c_i\in(h_i,2b_i)$ of $[b_i(c_i-b_i)]^3 c_i=729\cdot2^{4i}$. The root lies strictly above $h_i$ because $4b_i(h_i-b_i)<h_i^2$, and strictly below $2b_i$ since $h_i<2b_i$. For offsets at most $2^{4m}b_i$ use the general pair-sum bound. For larger offsets through $2^{4m}c_i$, their top binary unit is $2^{4m}b_i$ and the preceding product bound applies. The offset stays below $2^k$ in all seven cases. For $i=5$ the improvement in interval length over $6^{6/7}2^{4k/7}$ exceeds 19 percent.

## Four predecessor primes and cofactors

If $\omega(M)\le4$, one row takes at least three distinct predecessor primes. Indeed exactly one of $S\ge V,T\ge U$ holds by comparing $ST=N-1$ with $UV=N$. On that side all three bins are nonempty, since otherwise a positive normalized part would be at least its entire row sum. Consequently the other side takes at most one predecessor prime.

If $A$ takes at most one, at least one of $x,y$ is a multiple of $US$, hence $h>US$. Also $U\le h/2$. From $27T^2\le U^3$ we get $M^2<h^3/54$, implying the stated bound. If $C$ takes at most one, two of its normalized parts are positive multiples of $T$, so $2T<U$. Now $(US)^2=xy(h-1)/(abc)\le h^2(h-1)/4$. Thus $16M^2<h^2(h-1)$. With at most three primes the small side takes none. If $T=1$, the last identity gives $h^2(h-1)\ge4U^2 M^2\ge36 M^2$. If $S=1$ then $h\ge2U$, $h^3\ge216M^2$, which suffices as $h\ge2$.

At even weight, when $A$ takes at most one, $US\le h-1$ and the same argument gives $54M^2\le(h-1)^2h$. When $C$ takes at most one, either $A$ takes 2 and $ab$ is even, or $C$ takes 2 and its two normalized multiples of $T$ have different 2-adic valuations, giving $U>3T$. This gives the constant 32. With at most three primes, if $T=1$ and $A$ takes 2 then $U\ge4$ and $ab\ge2$, giving 128. If $T=1$ and $C$ takes 2, its three normalized entries have distinct valuations, giving $U\ge7$. If $S=1$, we have $\ell\ge2$ from the parity argument and $h^3\ge216\ell M^2$.

In the two-leading-digit setting, if $\omega(M)\le4$, the full conjecture follows already from $Bz(B+z-1)\le4M^2$. If $C$ takes at most one, this is the bound $4M^2<xy(h-1)$ above. If $A$ takes at most one, use $US\le h-1$, $U\le\min(x,y)\le z$ to get $27M^2\le (B+z-1)^2z$, which also forces $4M^2<Bz(B+z-1)$.

There is also a quantitative version for **arbitrary** row pairs. Suppose $M=Rm$ with $\omega(m)\le4$, allowing $R$ and $m$ to share primes. Among the five pairwise coprime bins $d,e,t_1,t_2,t_3$, at least one divides $R$, even under the partial assignment. If $d$ or $e$ divides $R$, then $US<Rh$ because $x$ is a multiple of $US/d$ and $y$ of $US/e$. The companion bound gives $M^2<J^2D^3 R^2 h^3/54$. If $t_i\mid R$, then $T/R<DU$ by the corresponding companion part, giving $M^2<J^2D^2 R^2xy(h-1)$. In either case
$$D^3 J^2 R^2 h^2(h-1)>4M^2. \tag{2}$$
In particular the common endpoint divisor grows whenever $R h^{3/2}=o(M)$.

## The shape closest to the bound

Put $\rho=6^6M^4/[h^6(h-1)]$. If $\rho>1/2$, (1) forces $P=1$, so $w=a=b=c=\ell=1$. The strict inequality $N-h>0$ and $N=1+deT$, $h=1+de$ give $T=1+jU$ for a positive integer $j$: the third actual part $de(T-1)$ is divisible by $U$, coprime to $de$. The equations become
$$
U=t_1t_2+t_1t_3+t_2t_3,\quad T=t_1t_2t_3,\quad
(d-U)(e-U)=U^2-1,\quad h=U(d+e)=de+1.
$$
The parity restrictions above force $N,h$ odd and $2\mid S$. Indeed if $2\mid T$, then $\ell$ is even by the two companion valuations. Thus all $t_i$ are odd, $d,e$ have different parities, and $j$ is even. Every $t_i>1$ since $T>U$. Both $d,e>U$, since $(d-U)e=Ud-1>0$ and symmetrically.

**At least six distinct primes divide $M$ in this shape.** If each $t_i$ were a power of a single prime, put $\alpha_{ik}=[C_k/t_i]_{t_i}/t_i$ for $k\ne i$, using least nonnegative residues. Since $C_i=VT/t_i=1\bmod t_i$, the integer numerator in $\alpha_{ik}$ is the inverse of $t_k$ modulo $t_i$. Hence $\alpha_{ik}+\alpha_{ki}=1+1/(t_it_k)$ by the two inverse congruences. But at base equal to the prime below $t_i$, addition has no carry, even in the second block of $\log_p(t_i)$ digits. Thus the two off-diagonal entries in row $i$ sum to less than one. Summing over the three indices contradicts the three inverse-pair identities. At least one $t_i$ therefore has two prime factors. The two bins $d,e$ supply two more primes. This extends the [second-block reciprocity test](crt_reciprocity.md) to a row carrying only a portion of the predecessor factors.

For $8\mid N$ and $A$ taking 2 we have $w,c$ odd, $ab,\ell$ even, $8\mid U$. The only way to have $P<96$ is $w=c=1$, $ab=\ell=2$. Again $T=1+jU$ since $h-1=de$ and $N-h=de(T-1)>0$. After a permutation $U=2t_2t_3+t_1t_3+t_1t_2$ with all $t_i$ odd. The [determinant-two lattice argument](extremal_parts.md#congruences-at-the-next-two-thresholds) applies exactly to this equation and $T-jU=1$, giving $U\equiv2,4\bmod8$. This proves the factor 96.

## An additive margin at every weight

It suffices to consider sequences along which $\rho\to1$, since any bounded upper estimate for $(h-H_M)/H_M^{3/4}$ forces that convergence. In the unit shape, both ratios on the right of (1) tend to one. Thus $x,y\sim h/2$ and $U^2=xy/(h-1)\sim h/4$. Put $y_i=T/t_i$. Stability of AM–GM on the compact simplex gives $y_i\sim U/3$, hence $t_i\sim\sqrt{U/3}$.

Here is the needed arithmetic separation of the $y_i$. The positive integer
$$K=(t_1+t_2+t_3)U-9T=\sum_{\rm cyc}t_i(t_j-t_k)^2$$
is nonzero because the three coprime bins exceed one. Each bin is odd, so $K=0\bmod4$ and $U=3\bmod4$ (expand the sum of the three pair products modulo four). Then $(K+9)/U=t_1+t_2+t_3-9j$ is positive and $3\bmod4$. Thus $K\ge3U-9$. The identity
$$U^3-27T^2=\frac U2\sum_{i<k}(y_i-y_k)^2+
3\sum_{\rm cyc}y_i(y_j-y_k)^2$$
has coefficients asymptotic to $3U/2$ when collected by each squared difference. The squared differences sum to $\sum t_i^2(t_j-t_k)^2\sim\sqrt{U/3}\,K$. It follows that $U^3-27T^2\ge(3\sqrt3/2-o(1))U^{5/2}$, as in the [individual-part estimate](extremal_parts.md#the-additive-improvement). Squaring the AM–GM ratio in (1) now gives $\liminf h^{1/4}(1-\rho)\ge3\sqrt6$. Since $1-\rho\le1-(H_M/h)^7\le7(h-H_M)/h$, this proves the asserted additive margin.

## Counting near the lower bound

For $X\le N\le 2X$ and $h\le Y\le 2X$ put $B=Y^7/X^4$. Equation (1) gives $w^2(ab)^3c\ell^2\ll B$; an empty parameter range contributes nothing. Set $n=Uw$. Then
$$n(ae+bd)=h=cde+1,\qquad wT-c=jU>0.$$
For the latter relation use $N-h=de(wT-c)$ and $\gcd(U,de)=1$. The inequality $ae+bd\ge2\sqrt{ab(h-1)/c}$ gives $U\ll\sqrt{cY/(ab)}/w$. The companion bound gives $j\le w\sqrt{U/\ell}$ and $\min t_i\le\sqrt U/\ell^{1/6}$. With $a,b,c,w,f_1,f_2,f_3$ fixed, choose a smallest bin $r=t_1$ and $j$. There are
$$O\big(\sqrt{cY/(ab)}\,\ell^{-2/3}\big)$$
choices, including the constant number of possible bin orderings. This counts positive integral choices within the two separate upper bounds using the uniform upper bound for $U$ above.

Put $\delta=wr-jf_1$. The remaining bins satisfy
$$
(\delta t_2-j f_2 r)(\delta t_3-j f_3 r)= j^2 f_2 f_3 r^2+\delta c .
$$
Both factors and $\delta$ are positive, by $wT=jU+c$ and the formula for $U$ in the companion row. Thus the divisor bound gives $O_\epsilon(X^\epsilon)$ choices of $(t_2,t_3)$, fixing $U$. A second factorization
$$(cd-na)(ce-nb)=n^2 ab-c$$
has positive factors since $(cd-na)e=nbd-1>0$ and $(ce-nb)d=nae-1>0$. Again there are $O_\epsilon(X^\epsilon)$ choices, and $M=wdeT$ is fixed. All integers being factored are bounded by fixed powers of $X$.

Sum the bound for $(j,r)$ over $c\ll B/[w^2(ab)^3\ell^2]$. The sum of $\sqrt c$ costs at most a constant times this upper bound to the power $3/2$. The remaining series in $w,a,b,f_1,f_2,f_3$ converge, since the respective denominator exponents are $3,5,5,11/3,11/3,11/3$. Absorbing the divisor bounds into $X^\epsilon$ gives
$$
\#\{X\le N\le 2X\text{ admitting a counterexample with }h\le Y\}
 \ll_\epsilon X^\epsilon\,Y^{11}/X^6 .
$$
Taking $Y=C(2X)^{4/7}$ and summing dyadic blocks proves the count in the theorem. In particular, for every large $Q=2^k$ the entire interval $Q<N\le Q+C Q^{4/7}$ contains at most $O_{C,\epsilon}(Q^{2/7+\epsilon})$ possible exceptions. The same argument gives density one in $Q<N\le Q+Q^\alpha$ for each $4/7\le\alpha<3/5$.

## Finite exclusion of the unit shape

At every $3\le N\le79\,840\,884\,316\,361\,452\,506\,499\,572$, any counterexample would satisfy $h^6(h-1)\ge2\cdot6^6M^4$. The same bound holds at every weight whenever $h\le2\,939\,398\,104\,466\,764$.

For $P=1$, order the companion bins $r\le s\le t$. Dividing $rst=1+j(rs+rt+st)$ by $rst$ gives $j<r\le3j+1$. For each even $j\le1000$ and odd $r$ in that range, factor
$$((r-j)s-jr)((r-j)t-jr)=(jr)^2+r-j.$$
There are 4,761 companion triples. The reciprocal-residue argument excludes 66 with all three bins prime powers. For all others, factor $U^2-1$ to recover $d,e$ and keep only coprime bins, opposite parities for $d,e$, and a unitary split of $N$. Among the 1,262,204 remaining assignments, the two parts $Ud,Ue$ always share a binary digit. Since 2 is assigned to their row, this excludes every one. A direct loop solving for $t$ from $r,s$ agrees on all 94 companion triples with $j\le50$.

Consequently the unit shape needs $j\ge1002$. From $U^3\ge27(1+jU)^2$ we have $U>27j^2$. Also $d+e\ge4U$: the positive factors $d-U,e-U$ sum to at least $2\sqrt{U^2-1}>2U-1$. Since $d,e$ have opposite parity, $h\ge U(4U+1)$. Substitute $U>27\cdot1002^2$ and $N=1+(h-1)(1+jU)$ to obtain the two finite bounds.

Run `python3 check_pair_sum_extrema.py` for 452,655 exact row-and-part comparisons through weight 40, including 73,389 with a partially retained power. The tests include (1), (2) and the parity restrictions. With python-flint installed, run `python3 check_pair_sum_reconstruction.py 1000` to certify the 46-digit example and reproduce the finite exclusion.
