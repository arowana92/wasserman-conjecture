# Digit sums of multiples exclude all exponents at once

Let $p$ be a prime, $s_p(x)$ the sum of the base-$p$ digits of $x$, and $\mu_p(d)=\min_{v\ge1}s_p(vd)$. An integer $t$ is called sturdy in base $p$ when $s_p(t)=\mu_p(t)$. The terminology and algorithms for this minimum have been studied independently of Wasserman's conjecture. [Clokie, Molina Lovett and Shallit](https://doi.org/10.1016/j.tcs.2022.05.029)

**Theorem.** Let $t\ge1$ be coprime to $p$ and have at most two distinct prime divisors. If $s_p(t)-\mu_p(t)\le2$, the full three-row conjecture holds at $1+p^a t$ for every $a\ge1$.

More generally, if $t=AB$ with $A,B>1$ powers of distinct primes, the same conclusion holds whenever
$$ 1+s_p(t)<\mu_p(A)+\mu_p(B)+\mu_p(t). \tag{1} $$

**Proof.** For at most one prime divisor of $t$, the predecessor of the weight has at most two prime divisors, a [previously proved case](all_weights_endpoint_reduction.md). In the other case set $W=p^a$. The [three-factor assignment](endpoint_factors.md#when-n-1-has-three-distinct-prime-factors) forces a row avoiding $p$ whose parts, in some order, are $Wx,Wy,Wz+1$. They satisfy $A\mid x$, $B\mid y$, and $t\mid Wz+1$. All are multiples of a common divisor $U\ge3$ of $N=1+Wt$, which is coprime to $pAB$. In particular $UA\mid x$, $UB\mid y$, and $Ut\mid Wz+1$. Since avoidance gives an addition with no base-$p$ carries,
$$
 1+s_p(t)=s_p(x)+s_p(y)+s_p(Wz+1)
 \ge\mu_p(UA)+\mu_p(UB)+\mu_p(Ut)
 \ge\mu_p(A)+\mu_p(B)+\mu_p(t).
$$
Neither $A$ nor $B$ divides a power of $p$, so each of its minima is at least 2. Both assertions follow.

## Prime-power partners of a power minus one

Suppose $M=p^h-1>1$ is a prime power. For **every** prime power $B$ coprime to $p$ with
$$ B<p^{2h}+p^h+1, $$
and for **every** $a\ge1$, the conjecture holds at $N=1+p^a M B$. In particular any Mersenne prime $2^h-1$ works with every odd prime-power partner in this range.

To prove the digit bound, every positive multiple of $p^L-1$ has digit sum at least $L(p-1)$. Split it into blocks of $L$ base-$p$ digits and add the blocks. This preserves its residue modulo $p^L-1$, keeps it positive, and does not increase digit sum. Repeat until it has at most $L$ digits. The only possible positive multiple left is $p^L-1$ itself. Thus $\mu_p(M)=h(p-1)$.

If $B$ is a power of the same prime as $M$, use the at-most-two-factor predecessor theorem. Otherwise $t=MB<p^{3h}-1$, so $p^{3h}-1-t$ is another positive multiple of $M$. Its digits are the complements of the $3h$ digits of $t$, including leading zeros. It follows that $s_p(t)\le 2h(p-1)$. Also $\mu_p(t)\ge\mu_p(M)=h(p-1)$ and $\mu_p(B)\ge2$. This proves (1).

For example all odd prime powers $B<2^{62}+2^{31}+1$ work at every $1+2^a(2^{31}-1)B$. The proof uses no factorization of these weights.

## Divisors, shifted positions and non-sturdy examples

If $C\mid p^L-1$, then
$$ \mu_p(C)\ \ge\ \frac{L(p-1)}{s_p((p^L-1)/C)}. $$
Indeed multiplying a multiple of $C$ by $(p^L-1)/C$ produces a multiple of $p^L-1$. Digit sums are submultiplicative by ordinary multiplication followed by carrying. Consequently the same lower bound applies to every multiple $t$ of $C$.

Take $Q=p^b$, $h\ge2$, $C=1+Q+\cdots+Q^{h-1}$ and $L=bh$. The quotient $Q-1$ has digit sum $b(p-1)$, so $\mu_p(C)=h$. Choose integers $j_i\ge0$ independently for $0\le i<h$, with $j_0=0$, and put
$$ t=\sum_{i=0}^{h-1}Q^{i+h j_i}. $$
Its terms occupy distinct positions. Since $Q^h\equiv1\pmod C$, it is a multiple of $C$ with digit sum exactly $h$, hence sturdy. Whenever $t$ has at most two prime divisors, every exponent at $1+p^a t$ works. Neither the lengths of the zero gaps nor the outer exponent enter the test.

For instance take $p=2$, $b=1$, $h=31$, moving the bit at 14 up by 31 places and the bit at 2 up by 62 places. The cofactor is
$$ t=(2^{31}-1)\,8589950981=18446779260229107707. $$
Both factors are prime. Its 31 bits already achieve the minimum for all its multiples.

The two-unit allowance in the theorem also applies to genuinely non-sturdy integers. For $t=131071\cdot131149=17189830579$ both factors are prime and $s_2(t)=19$. The first factor is $2^{17}-1$, and $s_2(1755t)=17$. This proves $\mu_2(t)=17$, so every $1+2^a t$ satisfies the conjecture.

Run the standard-library verifier for the stated arithmetic, a direct computation of the minima on small moduli, and independent factorial-valuation checks of the three-factor condition.

~~~sh
python3 check_sturdy_cofactors.py
~~~
