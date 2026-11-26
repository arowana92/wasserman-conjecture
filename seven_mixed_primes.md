# Two sparse odd predecessor primes

**Theorem.** Suppose $p,q$ are distinct odd prime divisors of $N-1$ and both $s_p(N),s_q(N)\le7$. Then every three positive trinomial coefficients of weight $N$ have a common prime divisor. The primes and all digit positions are unrestricted.

At most five units at either prime already suffice. If $N$ is even, both digit sums are even and the [two-prime theorem through six units](equal_predecessor_gaps.md#two-sparse-primes) applies. Thus it suffices that $N$ is odd and both sums are seven. By the [large-base theorem](seven_large_bases.md) any counterexample would have $p,q<257$.

An avoiding row at either base has two even paired parts with the same lowest binary bit, and therefore carries at 2. The endpoint coefficient $N(N-1)$ is also even. Consequently both sparse primes must be avoided by the **same** nonendpoint row, since the other one must avoid 2. Its even parts at $p$ are $1+p^u,p^i(1+p^d)$ by the [units-place theorem](seven_large_bases.md#locating-the-units-place-digit), and similarly at $q$. The first parts $1+p^u,1+q^v$ are different by unique factorization. Thus each has to equal the other base's unshifted pair. In particular each of the ordered prime pairs $(p,q)$ and $(q,p)$ must allow
$$1+p^u=q^j+q^k,\qquad u\ge1,\quad k>j\ge1. \tag{1}$$

## A bound on the exponent

For any distinct odd primes $p,q<257$, the least $t>0$ with $p^t=-1\bmod q$, if it exists, is half the order of $p$ modulo $q$. In (1) we have $u/t$ an odd integer, and lifting at $q$ gives
$$q^j=q^{v_q(p^t+1)+v_q(u/t)}\le K u,\qquad K=q^{v_q(p^t+1)}.$$
Testing $1\le t<q$ over the 2,756 ordered pairs gives $K\le28,561$; pairs without $t$ have no solution.

Put $\Gamma=q^k p^{-u}-1=-(q^j-1)/p^u$, which is nonzero. We have $k<6u$ since $q^k<p^u$ and $\log p/\log q<6$. Matveev's explicit bound for two positive rationals gives
$$\log|\Gamma|>-1.4\cdot30^5\,2^{4.5}(1+\log(6u))\log p\log q.$$
All logarithms in this section are natural. The two heights are below 6, so
$$u\log p<3\cdot10^{10}(1+\log(6u))+\log(28561u).$$
Consequently $u<M=10^{14}$: at $M$ the right side is below $M$, and its quotient by $u$ decreases thereafter, whereas $\log p>1$. [Matveev's bound in rational form, Languasco–Luca–Moree–Togbé, Theorem 2.1](https://archive.mpim-bonn.mpg.de/id/eprint/5025/1/mpim-preprint_2023-34.pdf)

For $\theta=\log p/\log q$, the verifier finds a reduced rational $a/b$ with $2M<b<2^{70}$ and $|b\theta-a|<1/(2M)$. Since $u<b$, the integer $ua-kb$ is nonzero, and therefore $|u\theta-k|>1/(2b)>2^{-71}$. On the other hand $q^k>p^u/2$ and
$$0<u\theta-k=\frac{\log(1+(q^j-1)/q^k)}{\log q}<\frac{2Ku}{3^u}.$$
Already at $u=80$ the last bound is smaller than $2^{-71}$, and it decreases. Thus $u<80$.

The approximants are verified by rational inequalities. For $1\le x\le2$ use $t=(x-1)/(x+1)$ and sum $2\sum_{i=0}^{99}t^{2i+1}/(2i+1)$ for $\log x$, bounding the tail above by $2t^{201}/(201(1-t^2))$. Rescaling the primes by powers of two gives rigorous dyadic enclosures for all logarithms. Continued fractions of the midpoint suggest the approximants; direct comparison with both interval endpoints certifies them. The largest denominator used is 261,977,540,510,658,979.

## The complete one-way list

The verifier evaluates $1+p^u$ for $1\le u<80$ and tests the two-power representation using exact division. Independently it generates the sums $q^j+q^k$ through the same maximum value and intersects the two lists. The only ordered pairs $(p,q)$ are
$$(11,3),(29,3),(29,5),(47,13),(83,3),(89,3),
(107,3),(131,11),(149,5),(181,13),(251,3).$$
Every one has $p>q$. Hence both directions of (1) cannot hold, proving the theorem. Run `python3 check_seven_mixed.py` for the prime, exponent and interval checks. [Verifier](check_seven_mixed.py)

## A sparse prime dividing the weight

Suppose instead that $p$ is odd, $p\mid N-1$, $s_p(N)=7$, and $q\mid N$ is a prime with $s_q(N)\le5$. The full conjecture holds if **any** of $s_q(N)\le3$, $q\equiv1\bmod4$ or $N\equiv1\bmod4$ holds.

Indeed $N,q$ are odd. The $p$-avoiding row has two even pairs, one containing 1, and an odd three-unit part divisible by $p$. A row avoiding $q$ needs at least one singleton $q$-power, since it partitions at most five units. In the $p$-avoiding row neither even part could be such a power, and the odd part is divisible by the different prime $p$. Thus the companion must avoid both $q$ and 2. It has exactly one odd part. Its two even parts need at least two $q$-units each, so there must be five units and the odd part must be $q^k$. The even parts have the form $q^i(1+q^d),q^j(1+q^e)$ with $d,e\ge0$. Their 2-adic valuations must differ. For an even gap the valuation is one, and for an odd gap it is $v_2(q+1)$. Consequently $q=3\bmod4$ and the gaps have opposite parity. A binary carry-free sum with an even part of valuation one and an odd part has $N=3\bmod4$. In fact $k$ must be even, since the odd part has zero bit at position one.

In the remaining configuration this companion cannot avoid any odd predecessor prime. Its pure $q$-power is coprime to $N-1$, so the two even parts would both have to vanish at such a prime. Their plus binomials have odd gcd one because their gaps have opposite parity, also when a gap is zero.
