# An individual common prime growing with the weight

Uniformly for $y\ge2$, the number of weights $X\le N<2X$ at which some three positive trinomials have **no common prime greater than $y$** is
$$ O\!\left(X\frac{\log^2(2y)}{\log X}\right). \tag{A} $$
The same bound holds if two positive trinomials have no common prime greater than $y$ dividing $N(N-1)$.

Let $y=y(X)\ge2$ be any function with $\log(2y)=o(\sqrt{\log X})$. Outside $o(X)$ weights $X\le N<2X$, **every three positive triples of weight $N$ have a common prime divisor greater than $y$**. For example one may take
$$ y=\exp\!\left(\frac{\sqrt{\log X}}{\log\log X}\right). $$
For three rows one can use the primes above $y$ from the following list. Take the greatest prime at most $N$ and the prime factors of $N,N-1,\ldots,N-5$. All rows may repeat. The assertions are asymptotic and leave the all-weight conjecture open.

## The endpoint gcd of one row, with a uniform error

Write $L=\log X$, $M=N-1$, and $E=NM$. On all but $o(X)$ weights in the interval, every positive trinomial $B$ has
$$ G_B=\gcd(E,K_B)> X e^{-\sqrt L}. \tag{1} $$
We need a uniform error beyond the $X^{o(1)}$ factor in the [two-congruence count](two_predecessor_congruences.md#recovering-the-row-from-a-product). Restrict first to weights with $\tau(N)\tau(M)\le L^{10}$, where $\tau$ is the number of divisors. The excluded count is $O(XL^{-7})$ by Cauchy–Schwarz and $\sum_{n\le x}\tau(n)^2\ll x(1+\log x)^3$. In the reconstruction there, the factorization $W=Vd_3$ is the only divisor-choice loss in a fixed box. Since $W\mid NM$, any $W$ representing a retained weight has at most $L^{10}$ divisors. The CRT count in the other direction has absolute constants, and there are $O((1+L)^5)$ boxes. Taking $Q=2X$ gives $O(Y(1+L)^{15})$ retained weights with $G_B\le Y$, proving (1).

We also use $\sum_{n\le x}\tau(n)^j\ll_j x(1+\log x)^{2^j-1}$ for fixed positive integers $j$. Indeed $\tau(n)^j\le d_{2^j}(n)$, where $d_m$ counts ordered factorizations into $m$ positive factors. For a prime exponent $e$, encode $j$ choices of exponents from $0,\ldots,e$ by counting, for each subset of indices, the levels $1,\ldots,e$ attained by exactly that subset. Those $2^j$ counts sum to $e$ and recover the exponents. The summatory bound for $d_m$ follows by summing reciprocals of the first $m-1$ factors.

## If one row avoids a moderately small prime

Put $Z=\exp(\sqrt L/100)$. With $o(X)$ further exceptions, **simultaneously for every prime $p\le Z$**, every positive trinomial $A$ avoiding $p$ and every positive trinomial $B$ have
$$ \gcd(E,K_A,K_B)> e^{\sqrt L}. \tag{2} $$
Here $p$ can be any prime in the range, independently of its divisibility into $E$.

Let $m_p=\lceil\log_p(2X+1)\rceil$. Among all $p^{m_p}$ padded words the proportion with digit sum above $(1/2+\epsilon)(p-1)m_p$ or below $(1/2-\epsilon)(p-1)m_p$ is at most $e^{-2\epsilon^2m_p}$ for each specified tail, by the [bounded-digit exponential moment](several_carry_counts.md#a-criterion-allowing-different-primes-and-row-lengths). Since $m_p\ge100\sqrt L$ and $p^{m_p}\le3pX$, the number of weights with digit sum above $(3/5)(p-1)m_p$ for **some** $p\le Z$ is $O(X Z^2 e^{-2\sqrt L})$. Discard them.

Choose a largest part $a$ of $A$. The interior binomial $C=\binom Na$ divides $K_A$, and is carry-free in base $p$. At least one of its two parts $z\in\{a,N-a\}$ has digit sum at most $(3/10)(p-1)m_p$, and both are at least two for $N\ge4$. The union $\mathcal Z$ of possible $z$ over all these bases has size
$$ |\mathcal Z|\ll X Z^2 e^{-8\sqrt L}. \tag{3} $$

If (2) fails then $G=\gcd(G_B,C)\le e^{\sqrt L}$. Set $R=G_B/G=US$, where $U\mid N$ and $S\mid M$. By (1), $R>Xe^{-2\sqrt L}$. At each prime of $R$, its exponent is at most $v_q(E)-v_q(C)$. The [partial-power residue lemma](endpoint_pair_bound.md#retaining-partial-prime-powers) applied to the two parts therefore gives $U\mid z$ and $S\mid z(z-1)$. For each $z$ there are at most $F(z)=\tau(z)^2\tau(z-1)$ choices of $(U,S)$, each allowing at most $1+X/R\ll e^{2\sqrt L}$ weights in one residue class modulo $R$.

Two applications of Cauchy–Schwarz and the preceding fixed moments give
$$
 \sum_{z\in\mathcal Z,\,z\ge2}F(z)
 \le |\mathcal Z|^{1/2}\left(\sum_{2\le z\le2X}F(z)^2\right)^{1/2}
 \ll |\mathcal Z|^{1/2} X^{1/2}(1+L)^{O(1)}.
$$
Indeed the sum of squares is bounded by the geometric mean of the eighth and fourth divisor moments. By (3), even after multiplying by $O(e^{2\sqrt L})$ the count is $O(Xe^{-\sqrt L})$ for large $X$. This proves (2).

## Removing the carry condition on the pair

The count of weights at which $N(N-1)$ has no prime factor in $(y,Z]$ is
$$ \ll X\left(\frac{\log(2y)}{\sqrt L}\right)^2. \tag{4} $$
Assume $6\le y<Z$. Increasing a smaller $y$ to 6 costs an absolute factor, and for $y\ge Z$ the bound is trivial. For $y<p\le Z$ forbid the two residues $0,1\bmod p$. A set of $j$ distinct primes with product $d$ imposes $2^j$ disjoint residue classes, containing $2^j X/d+O(2^j)$ weights. Truncate inclusion-exclusion after the even level $m=2\lceil L^{1/4}\rceil$. This is an upper bound by the Bonferroni inequalities. The combined error from the integer counts is $O((2Z)^{m+1})=O(X^{1/2})$ for large $X$.

Put $H=\sum_{y<p\le Z}2/p$. Partial summation of the [elementary prime sum](two_predecessor_congruences.md#optimal-density-one-exponent-for-one-row) gives $H=2\log(\log Z/\log y)+O(1)$ uniformly. The difference between the truncated normalized main term and $\prod_{y<p\le Z}(1-2/p)$ in absolute value is at most $\sum_{j>m}H^j/j!$. This is $O(L^{-3})$, since $H=O(\log L)$ and $j!\ge(j/e)^j$. The product is at most $e^{-H}\ll(\log y/\log Z)^2$. This proves (4).

Outside the exceptions in (4), choose $p\mid N(N-1)$ with $y<p\le Z$ and consider an arbitrary pair $A,B$. If both coefficients are divisible by this $p$, it is already a common endpoint prime exceeding $y$. Otherwise one row avoids $p$ and (2) applies. Write $m_y(n)$ for the product of the full powers of primes at most $y$ dividing $n$. To bound the weights with $m_y(N)m_y(M)> e^{\sqrt L}$, use the [uniform smooth-part moment](growing_endpoint_divisors.md#from-a-common-divisor-to-a-common-prime) at threshold $e^{\sqrt L/2}$ on either integer. The exceptional fraction is $O(\exp[-\sqrt L/(4\log(2y))])=O(\log^2(2y)/L)$. Thus the divisor in (2) has a prime factor greater than $y$. All discarded sets depend only on the weight. The exceptional counts in (1), (2) and (4) now give (A) for pairs.

## Three unrestricted trinomials

Increase $y$ to at least 6. For any row having two parts $a,b\le4$, its coefficient $K$ satisfies $E\mid abK$ by cancelling the two factorials one step. Apply the unrestricted pair result to the other two rows. Their common endpoint prime exceeding $y$ also divides $K$.

Now suppose the second smallest part of every row is at least 5. For almost all weights the greatest prime $r\le N$ has $g=N-r\le N^{0.14}$. Apply Guth–Maynard's Corollary 1.4 at length $\lfloor X^{0.138}\rfloor$ and shift the interval to end at $N$; the total exception count is $O(X^{0.138}+X\exp[-cL^{1/4}])$. If $r$ is not a common prime, a row avoiding it has two smallest parts summing to at most $g$, since $r>N/2$. [Guth and Maynard](https://arxiv.org/pdf/2405.20552)

For each $0\le h\le5$, the smooth-part moment also gives $m_y(N-h)\le N^{1/25}$ outside $O(X\exp[-L/(50\log(2y))])$ weights per shift, by the moment at exponent $1/(2\log(2y))$. Suppose that the three rows collectively avoid all primes exceeding $y$ in the proposed list. Apply the [weighted prime-power inequality](weighted_obstruction.md#a-weighted-divisibility-lemma) at height 6, counting only primes greater than $y$. Removing the smooth parts on the left gives a lower bound of order $N^{21(1-1/25)}$ for the product of the three row contributions, as the six exponents $6,5,\ldots,1$ sum to 21.

Use the [height-six cutoffs](typical_certificates.md#when-all-second-parts-are-at-least-5). In the two ordinary rows give cutoffs $(2,2,2)$ if the least part is at least two, or $(1,2,3)$ if it is one, in increasing order of the parts. Each contributes at most $N^9$. In the row avoiding $r$ give cutoffs $(2,4,0)$ or $(1,5,0)$, respectively. Its contribution is at most $g^{15}$, since a part equal to one contributes exactly one. This gives at most $N^{18+15(0.14)}$, contradicting $21(1-1/25)>18+15(0.14)$ for sufficiently large $N$. Both new exceptional counts are $O(X\log^2(2y)/L)$, proving (A) for every triple simultaneously.

The verifier checks 114,827 exact coefficients, 93,272 assignments from the endpoint gcd ratio, and the two-residue sieve on a full period. Among the 155 weights from 6 through 160, the proposed list provides a common prime exceeding 11 for all three-row choices at 77 weights. Run the finite checks with the standard library.

~~~sh
python3 check_growing_common_primes.py
~~~
