# The weights $5t+1$ and $7t+1$

Let $\mathcal T$ consist of the products $t=p^a q^b$ of two prime powers, with $a,b\ge1$ and the primes allowed to coincide. A weight is exceptional if three positive trinomial rows of that weight have coefficients with gcd one.

1. **For every $t\in\mathcal T$, at least one of $5t+1$ and $7t+1$ satisfies the full conjecture.**
2. **For each $w=5,7$ separately, every interval of $H$ consecutive integers contains at most $O_\epsilon(H^\epsilon)$ exceptional weights of the form $wt+1$ with $t\in\mathcal T$, for every $\epsilon>0$.** The implied constant is independent of the position of the interval.

The first assertion uses residues modulo 3 and 4. The second uses all three necessary prime bases and the higher-dimensional slicing theorem of Corso and Shmerkin. The latter bound allows $t$ and both of its prime bases and exponents to vary.

## Incompatible residues

Suppose $5t+1$ is exceptional. The [factor-five proof](predecessor_five_seven.md#full-factor-five) gives $t=2^k u$, $k\ge3$, $u$ odd, $3\mid 5t+1$, and a binary carry-free cut $3t,2t+1$. In particular $t\equiv1\bmod3$. The bits of $3u$ and $2u$ are disjoint, so $u\equiv3\bmod4$. Indeed $2u$ has its bit at position 1 set, forcing $3u\equiv1\bmod4$.

For an exceptional odd weight $7v+1$ with $v\in\mathcal T$, the [factor assignments](carry_intersections.md#necessary-additions) give either $3\mid7v+1$ or the binary carry-free cut $5v,2v+1$. If $v$ is even with odd part $z$, the latter cut forces $5z\equiv1\bmod4$, or $z\equiv1\bmod4$. The former alternative requires $v\equiv2\bmod3$. Thus a factor-five exception at $t$ and a factor-seven exception at $v$ are incompatible whenever both inputs are even, $t\equiv v\bmod3$ and their odd parts agree modulo 4. In particular take $t=v$. For odd $t$ the factor-five theorem already suffices. The integer implication from the two low-layer inequalities and the mod-3 conditions is checked by Lean 4.12.0. [Formal residue proof](CarryOrder.lean)

An individual consequence is that $7v+1$ satisfies the conjecture whenever $v\in\mathcal T$ is even, its odd part is 3 modulo 4, and $3\nmid7v+1$.

## Three carry graphs

For the factor-five family the three required additions for one cut are those of $3t,2t+1$ in bases $2,3,5$. For the factor-seven family the even weights and the multiplier-six case among odd weights have monotone base-seven digits in $t$. Among the other odd-weight cases one of the last two rows of the following table applies. The two summands are $at,(w-a)t+1$. [Necessary additions](carry_intersections.md#necessary-additions) · [Monotone digits for multiplier six](predecessor_five_seven.md#a-restriction-on-every-base-seven-digit)

| $w$ | $a$ | Required bases | Upper bounds for the three dimensions | Sum |
|---:|---:|---|---|---:|
| 5 | 3 | 2, 3, 5 | 0.61, 0.69, 0.51 | 1.81 |
| 7 | 3 | 2, 3, 7 | 0.57, 0.66, 0.62 | 1.85 |
| 7 | 5 | 2, 5, 7 | 0.65, 0.74, 0.59 | 1.98 |

For each base $r$ take the [carry graph](carry_blocks.py) for the two summands, starting with carries $(0,1)$. A digit $v$ at state $(i,j)$ produces $av+i,(w-a)v+j$. Retain the transition if the residues modulo $r$ sum to less than $r$, then pass the quotients to the next state. Let $m=256$, $B_r=r^m$, and let $D_r$ be the values of all $m$-digit input blocks readable from any reachable state. The block value uses the usual place weights, with the first input digit least significant. Every full block of the scalar $t$ lies in $D_r$ if the cut is carry-free. The zero block belongs to the alphabet, since state $(0,0)$ is reachable by zeros.

The number of possible block digits is at most the total number $T_r$ of length-$m$ paths, counted from all starting states. The checker computes this twice, by successor sums and by powering the integer transition matrix. Each dimension bound $j/100$ in the table is certified by $T_r^{100}\le r^{mj}$. Let $A_r$ be the compact set of reals in $[0,1]$ whose fractional base-$B_r$ digits lie in $D_r$. Its dimension is at most $\log |D_r|/\log B_r$, by a cover by digit intervals. It is forward invariant under multiplication by $B_r$ modulo 1.

## From the digit sets to moving intervals

The three bases $B_r$ are pairwise multiplicatively independent by unique prime factorization. For the three sets $A_r$ in any row of the table the dimension sum is strictly less than 2. Corollary 1.17 of Corso and Shmerkin therefore covers an intersection of their affine images by $O_{\epsilon,K}(\delta^{-\epsilon})$ intervals of length comparable to $\delta$, uniformly in translations and in positive slopes bounded by a fixed $K$. [Corso and Shmerkin, *Dynamical self-similarity, $L^q$-dimensions and Furstenberg slicing in $\mathbb R^d$*, Theorem 1.15 and Corollary 1.17](https://arxiv.org/pdf/2409.04608)

Take any interval $[L,L+H)$ of positive integer inputs, $H\ge2$. For each base choose the smallest power $M_r$ of $B_r$ at least $H$, so $1\le M_r/H<B_r$. The quotient $b_r=\lfloor t/M_r\rfloor$ takes at most two values on the interval. The residue $(t-b_rM_r)/M_r$ belongs to $A_r$ because all its finite block digits are allowed and its expansion may be continued by zeros. For each fixed choice of quotients,
$$
 \frac{t-L}{H}\in\bigcap_r\left(\frac{M_r}{H}A_r+
                  \frac{b_rM_r-L}{H}\right).
$$
The slopes are bounded by the largest of the three fixed block bases. Distinct integers give points at least $1/H$ apart. Apply the uniform covering bound at $\delta=1/H$ and sum over the at most eight quotient choices. This proves $O_\epsilon(H^\epsilon)$ for each of the three cuts.

In the monotone base-seven case, choose $7^j\ge H$ minimally. The low $j$ digits of $t$ weakly decrease from low to high. There are at most $\binom{j+6}{6}$ such words, by counting the multiplicities of the seven digits. Distinct integers in the interval have distinct residues modulo $7^j$. This contributes $O((\log(2H))^6)$, which the stated bound absorbs. An interval of $H$ consecutive weights gives an interval of at most $H$ consecutive inputs for each $w$.

Run `python3 check_paired_power_orders.py` with python-flint to reproduce the dimension bounds. Run `lean CarryOrder.lean` with Lean 4.12.0 to check the finite residue implication.
