# Three primes for two rows at products of four nearby powers

Let $P,P+a,P+b,P+c$ be powers of four distinct primes, with $0<a<b<c$. Put $s=b-a$, $t=c-b$ and $F=P(P+a)(P+b)(P+c)$. Suppose
$$s(s+t)>a(2s+t),\qquad bs>t(b+s)\ \text{ or }\ t>b+s. \tag{1}$$
For **every** positive integer $d$ with $100dc^3\le P$, any two positive trinomial coefficients of either weight $dF$ or $dF+1$ share a prime from the bases of $P,P+a,P+b$. The full three-row conjecture follows at both weights by the [endpoint reduction](all_weights_endpoint_reduction.md). No coprimality condition is imposed on $d$.

For example $P=1{,}389{,}209$ and $(a,b,c)=(2,8,24)$ give two consecutive 25-digit weights. The same set of four primes certifies every three-row choice at both weights and has minimum possible size at each. A count below gives infinitely many products satisfying (1) without assuming prime tuples of any fixed shape.

## Carry freedom at two powers

Write $n=dF$ and choose any pair $W<Z$ of the first three powers. Set $u=Z-W$. We prove that a positive binomial row at $n$ cannot avoid both base primes. Carry freedom at a prime implies carry freedom in blocks of any power of that prime. Hence each part of a simultaneously avoiding row is divisible by $WZ$ and its quotient by $W$ has digits in base $W$ bounded by those of $n/W$.

At $W=P$ the digits of $n/W$ are the coefficients of $dX(X+a)(X+b)(X+c)/X$ at $X=W$. At $W=P+a$ the corresponding cubic is $d(X-a)(X+s)(X+s+t)$. Its constant is negative, its linear coefficient is strictly positive by (1), and its quadratic coefficient is positive because (1) implies $a<s$. Its base-$W$ normalization borrows one from the linear coefficient. Thus in both cases the digits in positions $1,2,3$ are bounded by $3dc^2,3dc,d$, respectively. The constant digit may be as large as $W-1$.

Write the quotients of the two parts by $W$ as $C_i(W)$, with their digit polynomials $C_i(X)$. Since $Z$ divides these quotients, $C_i(-u)$ is a multiple of $Z$. The digit bounds give
$$-7dc^3\le C_i(-u)<W+7dc^3.$$
Only $0$ and $Z$ are possible. In the latter case replace $C_i(X)$ by $C_i(X)+X-W$. The polynomial still takes the same positive value at $W$ and now vanishes at $-u$. Denote the result in either case by $(X+u)(\alpha X^2+\beta X+\gamma)$. Comparing coefficients down from degree three gives
$$0\le\alpha\le d,\qquad |\beta|\le4dc,\qquad |\gamma|\le8dc^2.$$

The quotient of the original part by $Z$ now equals $G_i(Z)$, where
$$G_i(Y)=(Y-u)\big(\alpha(Y-u)^2+\beta(Y-u)+\gamma\big).$$
It is a nonzero polynomial with a positive root $u$. It has a positive coefficient because $G_i(Z)>0$, so it must also have a strictly negative coefficient. Every coefficient has absolute value at most $19dc^3<Z/4-1$, by expansion.

At $Z=P+a$ the normalized total polynomial for $n/Z$ is the cubic with just one negative coefficient considered above. At $Z=P+b$ it is
$$d(Y-b)(Y-s)(Y+t)=d\big(Y^3+(t-b-s)Y^2+(bs-t(b+s))Y+bst\big).$$
In either case at most one of its four coefficients is nonpositive, and their absolute values are below $Z/4-1$. When normalizing signed coefficients in base $Z$, the incoming borrow at any place is zero or one. A strictly negative coefficient thus yields a digit greater than $Z/2$, while a strictly positive coefficient of the total yields a digit less than $Z/2$. A carry-free addition cannot accommodate negative coefficients from both parts at the same position, and cannot accommodate any opposite a strictly positive total coefficient. Two negative positions are needed and only one is available. This proves the binomial incompatibility for all three pairs.

At $n$ one can group two parts of any trinomial to get a binomial row. At $n+1$ a row avoiding both primes has residues $(0,0,1)$ up to permutation at each power. Some part is divisible by both. Group the other two and subtract one from their sum, which preserves positivity and carry freedom at both primes. Thus the binomial incompatibility at $n$ transfers to trinomials at both weights. Each trinomial row can avoid at most one of the three primes, so any two rows share one. [Transfer proof](clustered_predecessors.md)

## How often the gap conditions hold

Fix $0<\epsilon<1/4$. For large $X$ and integers $1\le D\le X^{1/4-\epsilon}$, there are
$$\gg_\epsilon \frac{X^{1/2}}{D(\log X)^4}$$
distinct products $F\le X$ for which the same three primes work at every $dF,dF+1$ with $1\le d\le D$. In particular at fixed $D=1$ the gap exponent $1/3$ applies to $\gg X^{1/2}/(\log X)^4$ pairs of consecutive weights.

Take $Y=X^{1/4}/2$ and partition $[Y,2Y]$ into intervals of length at most $H=\lfloor (Y/(100D))^{1/3}\rfloor$, using $O(Y/H)$ intervals. Here $\log H\gg_\epsilon\log Y$. By the prime number theorem the intervals containing $k\ge c_0H/\log Y$ primes account for $\gg Y/\log Y$ primes if the absolute constant $c_0>0$ is small enough. In such an interval, for any fixed $\eta>0$ there is a fixed $\lambda>0$, depending only on $\eta,\epsilon,c_0$, such that any interval of length $\lambda H$ contains at most $\eta k$ primes for sufficiently large $Y$. This follows from the upper bound $2L/\log L$ for intervals of length $L$. [Montgomery and Vaughan, equation (1.12)](https://personal.science.psu.edu/rcv4/personal/Publications/large_sieve.pdf)

Choose $\eta=1/8$. An upper-quarter prime is at distance greater than $\lambda H$ from every lower-half prime, by the primes between their ranks. Partition the original interval into $O(1/\lambda)$ cells of length at most $\lambda H/4$. One cell contains $m\gg_\epsilon k$ of the lower-half primes. Inside this set choose the third prime from its upper quarter. All its lower-half primes are at distance greater than $\lambda' H$ from this choice, by applying the same short-interval bound with a smaller fixed $\lambda'>0$ so that such intervals contain at most $m/8$ primes. Take the first two primes in order from one common cell of length at most $\lambda'H/3$ containing $\gg_\epsilon k$ of the latter lower-half primes. Finally take the fourth prime from the original upper quarter. Each such choice gives $s>2a$ and $t>b+s$, since the first three primes lie in an interval of length $\lambda H/4$. These inequalities imply (1). There are $\gg_\epsilon k^4$ quadruples per counted interval. Summing and using its lower bound on $k$ gives $\gg_\epsilon YH^3/(\log Y)^4$ distinct products by unique factorization. Each satisfies $100Dc^3\le Y\le P$.

## A shared minimum four-prime certificate

With the example above take $n=3{,}724{,}611{,}517{,}745{,}287{,}445{,}819{,}539$ and $r=n-2$. The certificate is
$$\{1{,}389{,}209,\ 1{,}389{,}211,\ 1{,}389{,}217,\ r\}.$$
More generally, any prime $r$ with $n+1-P\le r\le n$ and $2r>n+1$ completes the three primes to a shared certificate. A row avoiding $r$ has a part at least $r$, so its two remaining positive parts sum to at most $P$. Avoiding one of the first three primes would make at least one of these two parts divisible by its power $W\ge P$ at either weight. Thus each row can avoid at most one of the four primes.

A prime $q$ divides every trinomial at a weight $N$ exactly when the base-$q$ digit sum is at most two, since otherwise its digit units can be distributed among three positive parts. If $q\nmid N$ and $N>2$, such a digit sum would require $N-1$ to be a power of $q$. Here both $n-1$ and $n$ have at least two distinct prime factors. For $q\mid N$ all digit sums exceed two, checked using the four factors of $n$ and the exact decomposition
$$n+1=2^2\cdot5\cdot42683\cdot4363108869743560019.$$
Therefore every prime can be avoided by some row. Any proposed certificate of at most three primes is defeated by taking one avoiding row per listed prime, so four is minimum at each weight.

Run `python3 check_quartic_triangle.py` with python-flint installed. It uses recursive order certificates for all required primes. Independently of the polynomial proof, the code enumerates all higher digits allowed at $W$. Divisibility by $Z>W$ then determines the units digit uniquely if it exists. It tests the ensuing binomial cuts by exact factorial-floor valuations at both primes. The three pairs give 13,496, 2,896 and 2,990 proper cuts to test, respectively, with zero simultaneous avoidances. On 238 smaller instances the same enumeration agrees with direct integer binomial coefficients, including 71 simultaneous avoidances.
