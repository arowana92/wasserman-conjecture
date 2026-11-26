# Counting a trinomial from two predecessor congruences

For $N\ge3$ put $M=N-1$ and let $H(N)$ be the minimum of $\gcd(NM,K_B)$ over positive trinomial rows $B$ of weight $N$. Uniformly for $1\le Y\le Q$, on $Q/2\le N<Q$,
$$ \#\{N:H(N)\le Y\} \ \ll\ Y Q^{o(1)}. \tag{1} $$
In particular, for each fixed $\delta>0$, every trinomial coefficient has a divisor in common with $NM$ exceeding $N^{1-\delta}$, outside $O(X^{1-\delta+o(1)})$ weights $N\le X$ for $0<\delta<1$. Exponent one is optimal for a density-one assertion. For each $\delta>0$ a positive lower proportion of weights have $H(N)<N^{1+\delta}$.

The number of possible exceptions to Wasserman's full three-row conjecture through $X$ is **$O(X^{0.8689})$**. The fraction of possible exceptions is $O(X^{-1/10000000})$ uniformly in intervals of length at least $X^{29/40}$ in $[X,2X]$. For a binomial and a trinomial (the binomial parts must both be at least two), all such pairs share an endpoint prime outside $O(X^{0.9375})$ weights.

## Recovering the row from a product

Use any divisor $V$ of $N$ dividing all three parts of $B$, and any divisor $T$ of $M$ such that, for each prime power exactly dividing $T$, the three actual part residues are $0,0,1$ in some order. Put $U=N/V$, $S=M/T$, $K=US$. We bound the number of weights admitting **any** such data with $K\le Y$ by the right side of (1). There is no coprimality requirement within either factorization $UV=N$, $ST=M$. The [partial-power lemma](endpoint_pair_bound.md#retaining-partial-prime-powers) lets us take the removed powers to be exactly $\gcd(NM,K_B)$, giving (1).

Group the powers of $T$ according to which part has residue one. Write the resulting pairwise coprime products in increasing order $d_1\le d_2\le d_3$, so $T=d_1d_2d_3$. Empty groups give factors 1. The positive normalized parts are $z_i T/d_i$ with
$$
 1\le z_i\le U d_i/T,\qquad V z_i T/d_i\equiv1\pmod{d_i}. \tag{2}
$$
Put all sizes $U,S,d_1,d_2,d_3$ in dyadic boxes. In the estimates below size letters on right sides may denote the box scales, up to absolute factors. Write $D=d_1d_2\le T^{2/3}$.

Fix $d_1,d_2,z_1,z_2$. The single integer $W=Vd_3$ then obeys
$$
 Wz_1 d_2\equiv1\pmod{d_1},\qquad Wz_2 d_1\equiv1\pmod{d_2}. \tag{3}
$$
These fix one class modulo $D$, with vacuous congruences modulo 1. There are $O(1+Vd_3/D)$ possible $W$. Factoring $W$ gives $V,d_3$ up to $Q^{o(1)}$ choices. The third congruence in (2) fixes $z_3$ modulo $d_3$. Since $U/T=K/M=O(1)$, this allows $O(1)$ positive $z_3$. All three parts and $N$ are now known. Multiplying the counts of the four fixed parameters and the remaining choices bounds the count in the box by
$$
 Q^{o(1)}(U/T)^2 D^2(1+Vd_3/D)
 \ll Q^{o(1)}\big(Q^{-2/3}K^2 S^{-4/3}+K\big). \tag{4}
$$
Indeed $D^2\le T^{4/3}$ and $U^2V/T=NK/M$.

For completeness a second count fixes $S,V,d_1,z_1$. With $t=T/d_1$, the congruence $Vz_1t\equiv1\bmod d_1$ implies that $N=1+Sd_1t$ is fixed modulo $Sd_1^2$, and it is zero modulo $V$. These moduli are coprime. Thus the box contains at most
$$
 O\big(SVd_1(Ud_1/T)(1+Q/(SVd_1^2))\big)
   \ll Q^{2/3}S^{4/3}+K \tag{5}
$$
weights since $d_1\le T^{1/3}$. The geometric mean of the first terms of (4) and (5), omitting $Q^{o(1)}$, is $K$. The smaller of the two bounds is therefore $O(KQ^{o(1)})$. Sum over the logarithmically many boxes in each coordinate.

## Optimal density-one exponent for one row

It suffices to take $0<\delta<1/2$. Put $a=1-\delta/2$, $b=1-\delta/4$. For $X^a<p\le X^b$ prime take $N=1+mp\in[X,2X)$; then $2\le m<p$ for sufficiently large $X$. The row $(1,p,(m-1)p)$ has coefficient $N\binom{mp}{p}$. The binomial factor is prime to $p$ by its base-$p$ digits. Hence the endpoint gcd is at most $Nm<N^{1+\delta}$. Different such primes produce disjoint sets of weights since $X^{2a}>2X$. Each gives $\gg X/p$ choices of $m$.

Here $\sum_{X^a<p\le X^b}1/p$ is bounded below by a positive constant. One elementary justification uses $\sum_{p\le x}(\log p)/p=\log x+O(1)$. Indeed if $\psi(x)$ sums $\log p$ over prime powers $p^j\le x$, then $\psi(2m)-\psi(m)\le\log\binom{2m}{m}\le2m\log2$ by factorial valuations. Summing at powers of two gives $\psi(x)=O(x)$. Expand $\log(\lfloor x\rfloor!)$ as $\sum_{p^j\le x}\lfloor x/p^j\rfloor\log p$. Replacing the floors costs $O(x)$, the terms with $j\ge2$ after division by $x$ are bounded by $\sum_{n\ge2}\log(n)/(n(n-1))$, and the integral estimate for the factorial gives the asserted formula. Subtract at $X^a,X^b$ and divide by $b\log X$.

## The unrestricted three-row count

By the [endpoint reduction](all_weights_endpoint_reduction.md), an exception supplies two trinomials with no common endpoint prime. Choose the odd coefficient as row $A$. The weight and predecessor powers assigned to the other row $B$ give exactly the hypotheses of (2), with $U,S$ assigned to $A$. The count for $US\le Q^k$ is now $Q^{k+o(1)}$. For $US>Q^k$ use the [two-ended binary entropy bound, equation (4)](common_factor_count.md#information-that-every-individual-part-supplies). The following rational parameters give the bounds stated above; the verifier checks affine majorants over the full unit square.

| Length of $A$ | Length of $B$ | $k$ | $\lambda$ | $r_0$ | $s_0$ | Upper exponent |
|---|---|---:|---:|---:|---:|---:|
| 3 | 3 | .86884 | .72186 | .64567 | .72802 | .8689 |
| 2 | 3 | .93743 | .62837 | .61750 | .64982 | .9375 |

For a mixed pair in the opposite orientation the [bound of exponent .919](common_factor_count.md) still applies. For fixed high and low binary strings leaving $\lfloor\theta\log_2 Q\rfloor$ consecutive bits free, take $\theta=29/40$, $k=\theta-3/10000$. Apply the [prescribed-bits bound](common_factor_count.md#uniformity-with-prescribed-bits-at-both-ends) with the free-word density cutoff $1/2+1/2500$. Writing $h$ for binary entropy, the typical exponent is at most
$$ (1-\theta)h(1/3)+\theta\big[h(1/6)+\log_2(5)/7500\big]+(2/3)(\theta-k). $$
This exponent, $k$ and the tail exponent $\theta(1-2/(2500^2\log 2))$ are strictly below $\theta-10^{-7}$. Covering intervals by aligned binary blocks gives the interval assertion.

## Larger common divisors for pairs and triples

For two positive trinomials, if one coefficient is odd, then almost always their common divisor with $NM$ exceeds $N^{3/8}$. If an odd coefficient comes from a row with every part at least two, this bound increases to $N^{1/2}$. For **arbitrary** pairs, almost always either the common divisor with $NM$ exceeds $N^{4/25}$ or both coefficients are divisible by the same power of two exceeding $N^{4/25}$.

There is also a three-row conclusion allowing one row with unrestricted carries. For almost all $N$, any three positive trinomials, two of which have every part at least two and each have $v_p(K)\le(\log_p N)/10000$ at a prime of their own choice in $\{2,3\}$, have a common divisor with $NM$ exceeding $N^{1/256}$. All four assertions here hold outside $O(X^{1-10^{-8}})$ weights. In the statements that give a gcd with $NM$, the [smooth-part estimate](growing_endpoint_divisors.md#from-a-common-divisor-to-a-common-prime) also gives a common prime larger than any prescribed $X^{o(1)}$ on a set of density one.

To prove these claims, designate an arbitrary unrestricted trinomial $B$ and put $G_B=\gcd(NM,K_B)$. By (1) discard weights with $G_B\le N^{1-\rho}$ for any such row, taking $\rho=1/2000$. If $G=\gcd(NM,K_B,K_1,\ldots,K_j)\le N^\gamma$, distribute the powers in $G_B/G$ to the $j$ remaining rows, each to a row of least truncated valuation at that prime. Each assigned power at $q$ has exponent at most $v_q(NM)-v_q(K_i)$ in its assigned row, so the same part congruences hold as in the [several-row entropy criterion](several_carry_counts.md#a-criterion-allowing-different-primes-and-row-lengths). Their product exceeds $N^{1-\rho-\gamma}$, and one row therefore receives a modulus at least $N^{b}$, $b=(1-\rho-\gamma)/j$. Apply equation (4) of that criterion to this row. For unrestricted parts in a row with a specified valuation use its largest part as an interior binomial cut; divisibility transfers the valuation hypothesis and any common divisor to the original coefficient. For the proof with cuts use $G$ defined from their binomial coefficients. A failure for the original rows implies the same upper bound on this $G$.

The parameters below give a strict saving in the digit-entropy exponent and the upper digit-sum tail. In the third row of the table, failure of the power-of-two alternative supplies the row with the indicated valuation bound. In all four rows use $\epsilon=1/1000$ except the last, where $\epsilon=1/10000$.

| $j$ | length used | allowed $p$ | upper $v_p(K)/\log_p N$ | $\gamma$ |
|---:|---:|---|---:|---:|
| 1 | 2 | 2 | 0 | $3/8$ |
| 1 | 3 | 2 | 0 | $1/2$ |
| 1 | 2 | 2 | $4/25$ | $4/25$ |
| 2 | 3 | 2, 3 | $1/10000$ | $1/256$ |

The verifier checks 48,590 exact trinomial coefficients, 4,364 constructed three-bin rows and the [explicit row at $N=222{,}788{,}354$](endpoint_factors.md#an-obstruction-to-bounding-a-single-rows-avoided-product), for which the actual endpoint gcd is less than $N$. Run the reconstruction and exponent checks with the standard library.

~~~sh
python3 check_two_predecessor_congruences.py
~~~
