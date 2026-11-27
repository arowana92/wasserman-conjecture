# A common prime of the weight and a second one of its predecessor

Outside $O(X^{1-10^{-8}})$ weights $N\le X$, the following holds simultaneously for all positive trinomial rows. Let $A_1,A_2$ have coefficients $K_1,K_2$. For each row separately choose $p_i\in\{2,3\}$ such that
$$
 v_{p_i}(K_i)\le 10^{-4}\log_{p_i}N.
$$
Then $\gcd(N,K_1,K_2)>N^\gamma$. If all six parts are at least two, also $\gcd(N-1,K_1,K_2)>N^\gamma$. The exponents are

| Chosen bases | $\gamma$ |
|---|---:|
| 2 and 2 | .0490 |
| 2 and 3 | .0269 |
| 3 and 3 | .0048 |

In particular, with all parts at least two, there are **two distinct common primes** of $K_1,K_2$ and the endpoint coefficient $N(N-1)$, one from each consecutive integer. On a set of weights of density one, each of the two primes can exceed any prescribed $X^{o(1)}$ on $[X,2X]$. At weights divisible by six the restriction on parts holds automatically if each selected prime is completely avoided. Kummer's criterion then forces its prime to divide every part of the corresponding row.

## Retain one side at a time

Write $H=N$ or $N-1$ and $J=\gcd(H,K)$ for a single row with parts $Z_1,Z_2,Z_3\ge2$. With $p$ and the valuation bound as above, put $\alpha_2=.5245,\ \alpha_3=.5024$. We show that outside the stated exceptional set every such row satisfies $J>N^{\alpha_p}$ for **each** choice of $H$.

Let $S=H/J$. For $H=N$, $S$ divides every $Z_i$, since $K Z_i/N$ is an integer multinomial coefficient. For $H=N-1$, $S$ divides $Z_i(Z_i-1)$ for every $i$. In fact, at $q^e\parallel N-1$, if the part residues modulo $q^j$ sum to more than one, then their sum is more than one at every level from $j$ through $e$. Each such level contributes at least one to $v_q(K)$ by the factorial-floor formula. Thus the residues at $q^{e-\min(e,v_q(K))}$ are $0,0,1$ in some order whenever this power exceeds one.

On $X\le N<2X$ with a failure $J\le N^{\alpha_p}$, choose a witnessing row deterministically. For fixed $S$, the congruence $N\equiv0$ or $1\bmod S$ leaves $O(X^{\alpha_p})$ possible weights, because $S\gg X^{1-\alpha_p}$. Conversely each individual $Z_i$ determines $S$ up to $X^{o(1)}$ divisor choices. Give the weights in any stratum equal probability and write $\mathcal H$ for Shannon entropy in natural units. Conditional subadditivity gives
$$
\mathcal H(N)\le\frac13\sum_i\mathcal H(Z_i)+\frac23\mathcal H(N\mid S)+o(\log X).
$$
Indeed $\mathcal H(Z_i)=\mathcal H(S)+\mathcal H(Z_i\mid S)-\mathcal H(S\mid Z_i)$ and $\sum_i\mathcal H(Z_i\mid S)\ge\mathcal H(N\mid S)$ since $N=\sum Z_i$. This adapts the endpoint entropy argument to either consecutive integer individually. [Entropy argument for several rows](several_carry_counts.md)

Let $L=\lceil\log_p(2X+1)\rceil$, $\epsilon=10^{-4}$, and $s_p$ denote digit sum. Discard the weights with $s_p(N)>(p-1)(1/2+\epsilon)L$. The bounded-digit exponential moment bound gives $O_p(X^{1-2\epsilon^2/\log p})$ such weights by counting all padded words. For the rest,
$$
\sum_i s_p(Z_i)=s_p(N)+(p-1)v_p(K)\le(p-1)(1/2+2\epsilon)L.
$$
For $0<u<1$, the entropy of a digit $d$ is at most $\log(\sum_{j=0}^{p-1}u^j)+\mathbb E[d]\log(1/u)$ by relative entropy. Sum over positions and parts. Writing $\mu=(p-1)(1/2+2\epsilon)/3$, the exponent bounding the number of remaining exceptional weights is at most
$$
\frac{\log(\sum_{j=0}^{p-1}u^j)+\mu\log(1/u)}{\log p}+\frac23\alpha_p+o(1).
$$
Take $u=1/5$ at 2 and $u=29/100$ at 3. Both bounds before the $o(1)$ are below $1-2\cdot10^{-8}$, and both digit-sum tails save more than $10^{-8}$. Summing dyadic intervals proves the single-row statement. A row with a part equal to one always has $N\mid K$, so the weight-side statement extends to those rows too.

For a pair, the two divisors $\gcd(H,K_i)$ have lcm at most $H\le N$. Their gcd therefore exceeds $N^{\alpha_{p_1}+\alpha_{p_2}-1}$, giving the table. For the individual primes, if all prime factors of this gcd were at most $y$, the part of $H$ supported on primes through $y$ would exceed $N^\gamma$. The smooth-part bound applies separately to $N$ and $N-1$, giving relative density tending to zero when $\log(2y)=o(\log X)$. [Smooth-part estimate](growing_endpoint_divisors.md#from-a-common-divisor-to-a-common-prime)

## Counting odd coefficients coprime to the predecessor

For weights $3\le N\le X$, the total number of ordered rows with an odd coefficient coprime to $M=N-1$ is at most
$$
X^{h(1/4)+o(1)},\qquad h(1/4)=0.811278\ldots,
$$
where $h$ is binary entropy. This count allows arbitrary factorizations of $M$. It complements the computational exclusion through two billion. [Predecessor census](odd_coefficient_predecessor.md)

Choose $L=\lceil\log_2(X+1)\rceil$. Binary carry freedom partitions the $L$ positions into the three part masks and the unused mask $2^L-1-N$. One mask has at most $L/4$ ones, so it has at most $2^{Lh(1/4)}$ possible values by the binomial sum bound. If it is the unused mask, it fixes $N$. Otherwise it gives a part $z\ge2$. A part equal to one would use the unique exceptional position at every full predecessor prime power, forcing the other positive parts to be multiples of $M$ even though they are smaller. Thus $M\mid z(z-1)$ leaves at most $X^{o(1)}$ possibilities for $N$. For each $N$, assigning the exceptional position at each full predecessor power fixes the three parts modulo $M$. Each is strictly between zero and $M$, so there are at most $3^{\omega(M)}=X^{o(1)}$ ordered rows. These estimates use the usual uniform divisor bound. The choice among the four masks contributes only a constant factor.

## Exact small examples

The two primes in the density theorem need not both exist at a particular weight. For interior rows with both coefficients odd, the first failure for the $N$ side is at 511, for rows $(7,56,448)$ and $(73,146,292)$. The first for the $N-1$ side is at 127, for $(7,8,112)$ and $(9,36,82)$. If each row may instead avoid either 2 or 3, the first failures are 216 and 53, respectively.

Run `python3 check_separate_endpoint_divisors.py`. It verifies the logarithm bounds by rational enclosures, computes the coefficients as exact binomial products, and checks the claimed first weights. The searches cover 31,080 odd interior rows through 511 and 9,640 interior rows avoiding 2 or 3 through 216. [Verifier and exact examples](check_separate_endpoint_divisors.py)
