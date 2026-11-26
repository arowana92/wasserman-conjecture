# A common endpoint divisor for three or four rows with few carries

Let $K_i$ be the multinomial coefficient of row $i$ of weight $N$. The following bounds hold for **every** choice of rows satisfying the indicated conditions, outside $O(X^{1-\eta})$ exceptional weights $3\le N\le X$. Each row can use its own prime $p_i$ from the indicated set. The valuation condition is $v_{p_i}(K_i)\le\alpha\log_{p_i}N$.

| Rows | Allowed primes | $\alpha$ | $\gamma$ | $\eta$ |
|---|---|---:|---:|---:|
| Three positive trinomials, parts unrestricted | 2, 3 | $1/1000$ | $1/16$ | $10^{-6}$ |
| Three trinomials, every part at least two | all primes through 23 | $1/1000$ | $1/64$ | $5\cdot10^{-7}$ |
| Four trinomials, every part at least two | 2, 3 | $1/10000$ | $1/128$ | $10^{-8}$ |

In each case the conclusion is
$$ \gcd\big(N(N-1),K_1,\ldots,K_k\big)>N^\gamma. \tag{1} $$
In particular, the first row of the table applies if each coefficient is indivisible by at least one of 2 and 3. On a set of weights of density one, all the indicated tuples also have an individual common prime divisor of $N(N-1)$ exceeding any prescribed $y=X^{o(1)}$ on $[X,2X]$. These are assertions about the coefficients of the chosen rows themselves and do not require one of the rows to be $(1,1,N-2)$.

## A criterion allowing different primes and row lengths

Fix $k\ge2$, lengths $t_i\ge2$, finite sets of allowed primes $\mathcal P_i$, parameters $\alpha_i\ge0$, and $0<\gamma<2$. In this criterion all parts are at least two. Use natural logarithms for entropy. Put
$$
 f_p(m;u)=\frac{\log(1+u+\cdots+u^{p-1})+m\log(1/u)}{\log p}\qquad (0<u\le1).
$$
Choose positive $\epsilon$ and real numbers $b_i\ge0$ with $\sum_i b_i\le 2-\gamma$. For every $i$ and every $p\in\mathcal P_i$ choose a real $0<u_{ip}\le1$. If
$$
 f_p\left(\frac{p-1}{t_i}\left(\frac12+\epsilon+\alpha_i\right);u_{ip}\right)
   +(1-1/t_i)\max(1-b_i,0)<1 \tag{2}
$$
for every such $i,p$, then (1) holds outside $O(X^{1-\eta})$ weights for every $\eta>0$ less than all the strict margins in (2) and all $2\epsilon^2/\log p$. The constants may depend on these fixed parameters.

Work on $X\le N<2X$ with $L_p=\lceil\log_p(2X+1)\rceil$. Discard weights with $s_p(N)>(p-1)(1/2+\epsilon)L_p$ for any allowed $p$, where $s_p$ is digit sum. Their number for each $p$ is $O_p(X^{1-2\epsilon^2/\log p})$. Indeed among all $p^{L_p}$ padded words the digits are independent and uniform. A single digit divided by $p-1$ and centered has mean zero and lies in $[-1/2,1/2]$. Convexity on this interval bounds its exponential moment at $s>0$ by $\cosh(s/2)\le e^{s^2/8}$. Take $s=4\epsilon$ in Markov's inequality.

For each remaining exceptional weight choose a deterministic witnessing tuple of rows and primes. Let $G$ be the gcd in (1). For every $q^e$ exactly dividing $N$ or $N-1$, put $\ell=\min(e,v_q(K_1),\ldots,v_q(K_k))$. Assign $q^{e-\ell}$ to a row of minimum valuation, breaking ties by index. Denote by $U_i,S_i$ the products assigned to row $i$ from $N,N-1$, respectively. The [partial-power lemma](endpoint_pair_bound.md#retaining-partial-prime-powers) applies to every row length. Indeed if at a level $j\le e$ the sum of part residues exceeds $h=0$ or $1$, the residue of $N$ on that side, it exceeds $h$ at every level through $e$, contributing at least $e-j+1$ to the valuation. Thus at each retained power the residues of the assigned row sum to $h$. In particular every part $z$ of that row satisfies
$$ U_i\mid z,\qquad S_i\mid z(z-1). \tag{3} $$
Since $G\le N^\gamma$, the product $\prod_i U_iS_i=N(N-1)/G$ is at least $X^{2-\gamma}/2$. For at least one $i$, $U_iS_i\ge 2^{-1/k}X^{b_i}$. Partition the weights by one such index and its chosen prime $p=p_i$. There are a bounded number of strata.

In a nonempty stratum, give the weights equal probability. Let $P=(U_i,S_i)$, $Z_j$ the $j$th part of row $i$, and $\mathcal H$ entropy. By (3), $Z_j\ge2$ specifies $P$ up to $X^{o(1)}$ choices by the divisor bound. Given $P$, the weight belongs to one residue class modulo $U_i S_i$, since $\gcd(U_i,S_i)=1$. Consequently
$$ \mathcal H(N\mid P)\le \max(1-b_i,0)\log X+O(1). $$
The row determines $N$ by summation, and $N$ determines $P$ by the choice of witness. As $\mathcal H(P\mid Z_j)=o(\log X)$ for every part, conditional subadditivity as in the [shared-label inequality](common_factor_count.md#information-that-every-individual-part-supplies) gives
$$
 \mathcal H(N)\le \frac1{t_i}\sum_j\mathcal H(Z_j)
        +(1-1/t_i)\mathcal H(N\mid P)+o(\log X). \tag{4}
$$

Legendre's formula gives $\sum_j s_p(Z_j)=s_p(N)+(p-1)v_p(K_i)$. Thus the sum of the $t_iL_p$ expected base-$p$ digits is at most $(p-1)(1/2+\epsilon+\alpha_i)L_p$. The entropy of a digit $D$ is at most
$$ \log(1+u+\cdots+u^{p-1})+\mathbb E[D]\log(1/u), $$
by nonnegativity of relative entropy against the distribution proportional to $u^d$, $0\le d<p$. Sum this bound over digits and parts, using subadditivity. The mean part entropy is at most $L_p\log(p)$ times the first term of (2). Now (4) bounds the logarithm of the stratum size, proving the criterion on dyadic intervals and then by summation.

## Parameters and parts equal to one

For the table take $\epsilon=\alpha$ and all $b_i=(2-\gamma)/k$. In the first case, for each trinomial of weight $N\ge4$ choose a largest part $a$. Then $2\le a\le N-2$ and $\binom Na$ divides the trinomial coefficient. Apply the criterion to the three binomial rows $(a,N-a)$, each of length two. This inherits the upper valuation bound and carries the common divisor back to the original three coefficients. The first result also applies directly to any three interior binomial coefficients with $2\le a\le N-2$. In the remaining cases use length three directly.

The verifier uses the trial values $u$ below in increasing order of the primes and encloses each logarithm by rational arithmetic.

| $k,t_i$ | $1000u$ | Upper bound in (2) |
|---|---|---:|
| 3, 2 | 335, 436 | .9990 |
| 3, 3 | 201, 291, 425, 518, 636, 676, 735, 757, 791 | .9988 |
| 4, 3 | 200, 290 | .9997 |

For the common prime assertion, if all prime factors of $G>N^\gamma$ are at most $y$, the product of powers of primes at most $y$ in at least one of $N,N-1$ exceeds $X^{\gamma/2}$. The [smooth-part moment bound](growing_endpoint_divisors.md#from-a-common-divisor-to-a-common-prime) adds at most $O(X\exp[-\gamma\log X/(4\log(2y))])$ exceptional weights on the interval. This bound is uniform over the tuples.

At individual weights endpoint divisors can differ from unrestricted common divisors. For instance, the three odd binomial coefficients at weight 15 with cuts 3, 5 and 7 have gcd 13, coprime to $15\cdot14$.

The finite search covers all 814,229 unordered positive rows of weights 3 through 1,500 whose coefficients are indivisible by 2 or by 3. At each such weight every three eligible rows have a common endpoint prime. The verifier also tests 501,239 multirow assignments against exact coefficients. Run all checks with the standard library.

~~~sh
python3 check_several_carry_counts.py
~~~
