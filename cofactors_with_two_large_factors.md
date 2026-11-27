# Growing prime bases when the cofactor has at most two large coprime factors

Write $L=\log H$, $W=\log\log H$ and take $1\le R\le W/10000$. Consider integers $n\in[H,2H]$ with **no three pairwise coprime divisors all exceeding $\exp(L/(100R))$**. With at most
$$ O\big(H\exp(-L/(12R))\big) $$
exceptions among these integers, Wasserman's full conjecture holds at **every** $1+np^e$ in the range
$$ p\text{ prime},\quad r=\max(1,\log p)\le R,\quad e\ge1,\qquad e\log p\le \frac{LW}{10000r}. \tag{1} $$
The constant and the asymptotic threshold are uniform in $R$. The statement includes cofactors divisible by the scale prime. Every product of two prime powers satisfies the divisor hypothesis. So does every $n$ whose quotient by its two largest full prime powers is at most $\exp(L/(100R))$, interpreting missing factors as 1.

For instance $R=W^{1/3}$ eventually gives all primes $p\le\exp(W^{1/3})$ and includes every scale $p^e\le\exp(LW^{2/3}/10000)$ in this base range. For any fixed finite set of prime bases of unrestricted size, taking $R$ fixed proves simultaneous validity at superpolynomial scales on a positive lower proportion of all cofactors. The algebraic constants are deliberately conservative.

## How the factorization controls the modulus

Use the [scaled endpoint assignment](growing_prime_power_steps.md) at $B=p^e$, with $Us=K$, $UV=Bn+1$ and $sT=n$. As in that proof, split $T$ into three pairwise coprime products $d_i$, placing each assigned power at the position where the complementary row has residue one. Its normalized parts $y_i$ are positive multiples of $T/d_i$ and sum to $U$. By the hypothesis on $n$, some $d_i\le\exp(L/(100R))$, so
$$ K=Us\ge sT/d_i\ge H\exp(-L/(100R)). \tag{2} $$
If $n$ has at most two distinct prime factors, one product is 1 and in fact $K\ge n$.

Assign to each exceptional $n$ a canonical witness for this scale. In the row avoiding $p$ the positive slopes $a_1,a_2,a_3$ sum to $n$ without a base-$p$ carry. Each actual part $z=Ba_i+\mathbf1_{i=3}$ is at least 3 and specifies $P=(U,s)$ up to $\tau(z)\tau(z(z-1))$ choices. Indeed $U\mid z$ and $s\mid z(z-1)$. Uniformly under (1), the logarithm of this multiplicity is at most $L/(100r)$ for sufficiently large $H$. Here are explicit margins. Put $M=B(2H+1)+1$. The elementary [divisor estimate](growing_prime_power_steps.md#counting-with-the-large-scale-removed) gives $\max_{j\le x}\log\tau(j)\le3\log x/\log\log x$ for all sufficiently large $x$. Apply it at $M,M^2$. Since $\log\log M\ge W$ and $\log M\le 2LW/(10000r)+O(1)$, the sum is at most $9\log M/W<L/(100r)$.

The pair $P$ fixes one class of $n$ modulo $K$ since $s\mid n$, $Bn\equiv-1\pmod U$, and $U$ is coprime to $Bs$. By (2) this permits at most $1+\exp(L/(100R))$ choices in the interval.

## A digit bound uniform in the prime

Use $l=\lceil\log_p(2H+1)\rceil$ digits with leading zeros, normalizing each digit by $p-1$. For a uniform normalized digit $X$, $0\le X\le1$ and $\mathbb E X=1/2$. Convexity on $[0,1]$ gives $\mathbb E\exp(\pm(X-1/2))\le\cosh(1/2)\le\exp(1/8)$, the last inequality following term by term from the exponential series. Among all $p^l$ possible inputs the fraction with digit sum greater than $3l/4$ is at most $\exp(-l/8)$ by the exponential moment bound for independent digits.

On the remaining exceptional inputs choose $n$ uniformly. The normalized digit sums of the slopes add to at most $3l/4$ by carry freedom. A digit with any law of mean $m$ has entropy at most $\log p-3/8+m$, by relative entropy against the law proportional to $\exp(-X)$ on the $p$ possible digits. Thus the average of the three slope entropies is at most $l(\log p-1/8)$, by subadditivity across positions. The [conditional entropy inequality](common_factor_count.md#information-that-every-individual-part-supplies) bounds the logarithm of the count by
$$
 l(\log p-1/8)+\tfrac23\log(1+\exp(L/(100R)))+L/(100r).
$$
Use $l\ge L/r$ and $l\log p\le L+O(1+r)$. Both this count and the count of atypical digit sums are $O(H\exp(-L/(10r)))$, since $1/8-1/150-1/100>1/10$ and $r^2=o(L)$ uniformly. There are at most $O(\exp(R)LW)$ pairs in (1). Its logarithm is $O(W)=o(L/R)$, so summing gives the claimed exponent $1/(12R)$.

## Size of the cofactor family

There are $\gg H/R$ integers satisfying the divisor hypothesis. For every $j\le\exp(L/(100R))$ take $n=j\ell$ with $H/j\le\ell\le2H/j$ prime. Each $j$ gives $\gg H/(j\log H)$ choices by the prime number theorem in dyadic intervals, uniformly since $H/j\ge H^{0.99}$. The prime $\ell>\sqrt{2H}$ is unique for large $H$. Of three pairwise coprime divisors of $n$, at least two divide $j$. Sum over $j$ using the harmonic series. For fixed $R$ this gives the positive proportion asserted above. Among products of two primes alone, the exceptional bound likewise has relative density zero, as already the cofactors $2\ell$ give $\gg H/\log H$ such integers.

The [Lean proof](ScaledCofactorBins.lean) verifies the modulus inequality from the three normalized part bounds. Run `lean ScaledCofactorBins.lean` with Lean 4.12 and `python3 check_cofactor_factor_bins.py` for 2,000 complete divisor comparisons, 77,525 assignments and 175 exact digit-sum tails.
