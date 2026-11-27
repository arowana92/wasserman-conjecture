# A unit added to one part of a scaled prime-power triple

Let $p$ be any prime and $M\ge1$. Take a triple of powers of $p$, absorb their common factor into $M$, and call their sum after normalization $s$. Thus at least one of the three powers is 1. Form $A$ by scaling this triple by $M$ and adding 1 to any one coordinate. Its weight is $N=sM+1$. For any positive trinomial row $B$ of weight $N$, put
$$G=\gcd(N(N-1),K_A,K_B),\qquad K_A=N!/\prod_i A_i!.$$

**Theorem.** Every such pair satisfies $G\ge3$. More precisely,
$$
\begin{split}
 G^6&\ge729s && (p\mid N),\\
 G^3&\ge27\big(s/p^{v_p(s)}\big)^2 && (p\nmid N).
\end{split}
$$
There is no restriction on the gaps between the powers, on $M$, or on the carries in $A$ and $B$. By the [endpoint theorem](all_weights_endpoint_reduction.md), a nonendpoint row of this form cannot occur in any counterexample to the full conjecture.

## The weight power remaining after carries

Write $v=v_p(N)$ when $p\mid N$ and let $\lambda=v_p(K_A)$. We show that the retained weight power $p^r$, $r=\max(v-\lambda,0)$, satisfies $p^{2r}\le s$ whenever $r>0$. The [retained-power lemma](endpoint_pair_bound.md#retaining-partial-prime-powers) shows that $p^r$ divides all three actual parts. Therefore $p\nmid M$, the shifted coordinate has power 1, and the other powers $p^i,p^j$ satisfy $r\le i\le j$. Indeed the shifted coordinate must be divisible by $p$, and the unshifted normalized powers are consequently positive multiples of $p$. We have $s=1+p^i+p^j$.

Put $t=\min(i,v)$. Because $N=1+sM$, we have $M\equiv-1\bmod p^t$. The floor-sum formula for $\lambda$ counts carries with multiplicity, one layer for each modulus $p^w$. At every layer $i<w\le v$ the contribution is at least one, since $Mp^i$ has nonzero residue there while the sum of the three residues is zero modulo $p^w$.

Now take $j<w\le i+t$. The residues of the two unshifted parts are $p^w-p^i$ and $p^w-p^j$. Their sum is at least $p^w$, and is strictly greater except at $p=2,\ i=j,\ w=j+1$. At $w>v$ this supplies a carry. At $w\le v$ the full residue sum is a multiple of $p^w$, so strict inequality supplies a *second* carry. This counts distinct contributions even when the two stretches of layers overlap. With $x_+=\max(x,0)$,
$$ \lambda\ge (v-i)_+ +(i+t-j)_+-\mathbf1_{\{p=2,\,i=j\}}.$$
For $p\ne2$ or $i<j$, this implies $r\le\min(i,j-i)$: subtract the first term from $v$ to get $t$, then subtract $(t-(j-i))_+$. Thus $p^{2r}\le p^j<s$. In the remaining case $r\le1$ and $s=1+2^{i+1}\ge5$ because $i\ge r>0$. This proves the claimed bound.

## The complementary row

Remove the common prime powers $D=\gcd(N,K_A,K_B)$ and $J=\gcd(N-1,K_A,K_B)$, so $G=DJ$. Assign the retained powers to the rows with the smaller valuations. Write $U,S$ for the products assigned to $A$ from $N,N-1$, respectively. The [product bound](projected_rows.md) for the other row gives
$$27(N-1)^2\le D^3J^2 U^3S^2. \tag{1}$$
Each prime power in $U$ divides all three parts of $A$. A prime other than $p$ dividing the unshifted entries would divide $M$, hence could not divide the shifted entry. Therefore $U$ is a power of $p$. If $p\mid N$ then $U$ is either 1 or the retained power just bounded. Every prime power in $S$ divides an unshifted entry since two entries have residue zero there. In this case $p\nmid S$ and $S\mid M$. Now $N-1=sM$ and (1) gives $D^3J^2\ge27s^2/U^3\ge27\sqrt{s}$.

If $p\nmid N$, then $U=1$. The part of $S$ prime to $p$ divides $M$, and its $p$-part divides $sM$. Consequently $S\mid Mp^{v_p(s)}$. Use (1) once again and $G^3\ge D^3J^2$ to get the second inequality. Both cases give $G\ge3$.

For the qualitative conclusion at $p\mid N$ there is a simpler direct argument. In a putative counterexample assign $p$ to $A$ if possible. Assigning it to $B$ would give $U=1$, although the positive normalized parts of $B$ sum to $U$. Thus $A$ adds without a carry at $p$ and receives the full weight power $P=p^v$. The two unshifted slots $Mp^i,Mp^j$ have runs of $v$ maximal digits starting at $i,j$ because $M\equiv-1\bmod P$. They cannot overlap. Hence $v\le i,j-i$ and $P^2\le p^j<s$, inconsistent with $27s^2\le U^3$.

## Binary offsets

In the [three-power diagonal theorem](three_power_diagonals.md), a shifted target gives this very row with $M=p^c+l$. The only remaining possibility is therefore
$$ C=l(1+p^h+p^k),\qquad l,\ lp^h,\ lp^k\text{ adding without a carry at }p.$$
In particular $s_p(C)=3s_p(l)$, so any offset whose base-$p$ digit sum is not a multiple of three excludes every high exponent triple in that theorem. The [sharper height bound and finite computations](three_power_offset_range.md) also cover small exponents at specific offsets.

For binary offsets $2\le C\le1,048,576$ two independent enumerations exclude all high exponents at 1,037,565 offsets. The remaining 11,010 offsets admit 13,810 offset/difference pairs. At most four pairs occur at one offset. One enumeration constructs the three shifted bit sets and the other factors each offset and checks its divisors.

~~~sh
python3 check_shifted_power_rows.py
~~~

The direct arithmetic checks compare endpoint gcds computed from exact binomial coefficients in 13,128 cases through weight 155, with duplicate companion gcds collapsed. The separate carry check covers 26,624 profiles, including 5,028 partially retained weight powers.

[Lean checks](ShiftedPowerCarries.lean) verify the deduction from the two layer counts to the retained exponent bound and the arithmetic at the reduced height threshold.
