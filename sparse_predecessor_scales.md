# The exact predecessor gcd along the direction (2,3,6)

Let $N=11V$ and let $K(V)$ be the coefficient of $(2V,3V,6V)$. At **every** positive multiplier,
$$
 \frac{N-1}{\gcd(N-1,K(V))}=
 \begin{cases}
 2 & V=(2^m+1)/3,\quad m\ge3\text{ odd},\\
 3 & V=(3^n+1)/2,\quad n\ge2,\\
 1 & \text{otherwise}.
 \end{cases}                                                         \tag{1}
$$
Thus the least uniform multiplier $E$ with $N-1\mid E K(V)$ is 6, but a single scale can lose at most a factor of 3. Since $11\mid K(V)$, its gcd with the endpoint is at least $11(N-1)/3$ at every weight on this ray.

## Sparse multiples forced by a lost power

Take a primitive direction with two entries $b,(p^h-1)b$ where $p$ is prime, $h\ge1$, and $p^k\parallel b$ with $k\ge1$. If the predecessor defect at $p$ is at least $j>0$, then
$$
              s_p((p^h-1)bV)\le(p-1)(h+k-j).                          \tag{2}
$$
Here $s_p$ is digit sum and $V$ must be a unit at $p$. In the [decrement formula](sharp_direction_predecessor.md) decrement the unit entry and $bV$. The multinomial on the three resulting parts contains the binomial on $bV-1,(p^h-1)bV$. This binomial equals $p^{-h}\binom{p^h bV}{bV}$, whose valuation is $s_p((p^h-1)bV)/(p-1)-h$ by Legendre's digit-sum formula. The decrement valuation is at most $k-j$, proving (2).

If losses in one fixed direction force (2) at two distinct primes, only finitely many scales realize both losses. More generally this applies to two fixed directions at common weights with respective losses at distinct primes. If the digit bounds are on $cV$ and $dW$, with $N=uV=vW$, apply the two-base digit-sum theorem to $cdN$. Indeed digit sums are submultiplicative, so at the first prime $s_p(cdN)\le s_p(du)s_p(cV)$ and likewise at the second. Stewart's effective lower bound on the combined digit complexity tends to infinity. [Stewart, Theorem 1](https://uwaterloo.ca/pure-mathematics/sites/default/files/uploads/documents/j-reine-ange-math-1980.pdf)

## Evaluating the two sequences

For $(2,3,6)$ the coarse pair-gcd bound is 6. At 2, equation (2) gives $s_2(3V)\le2$, with $V$ odd. The positive odd multiple $3V$ must be $2^m+1$ for an odd $m$. The binomial on $2V,6V$ has valuation $s_2(6V)=2$. For odd $m\ge5$ the binary ones of $V$ occur at zero and at $1,3,\ldots,m-2$. Thus $8V$ and $3V=2^m+1$ are disjoint in binary. The outer binomial has valuation zero, and $v_2(11V-1)=3$ since $3(11V-1)=11\cdot2^m+8$. At $m=3$ the outer binomial has valuation two and the predecessor valuation is five, again giving defect one. At $m=1$ there is no loss.

At 3 the required bound is $s_3(2V)\le2$. A positive even ternary integer prime to 3 with digit sum at most two is either 2 or $1+3^n$. The binomial on $3V,6V$ has valuation $s_3(6V)/2=1$ in the latter case. The base-three digits of $V=(3^n+1)/2$ are a final two preceded by $n-1$ ones. For $n\ge3$, adding $9V$ to $2V=3^n+1$ has no carry, and $v_3(11V-1)=2$. At $n=2$ the outer binomial has one carry and the predecessor valuation is three. The cases $V=1$ and $n=1$ give no loss.

The two nontrivial sequences cannot meet. Equality would give $2^{m+1}-3^{n+1}=1$. Since $m$ is odd, the two consecutive odd factors $2^{(m+1)/2}-1$ and $2^{(m+1)/2}+1$ would be coprime factors of a power of three. The smaller must equal one, forcing $m=1$. This proves (1) without any ineffective bound.

## Arbitrary losses and their frequency

For any $h\mid m$, $k\ge1$, set $b=p^k(p^m-1)/(p^h-1)$. On the primitive direction $(1,b,(p^h-1)b)$ the exact $p$-exponent of $E$ is
$$
                          \max(0,k+h-m).
$$
Every positive multiple of $p^m-1$ has base-$p$ digit sum at least $m(p-1)$. To see this, fold its digits into $m$ columns by exponents modulo $m$. Replace any $p$ units in a column by one in the next cyclic column. Each replacement preserves the residue modulo $p^m-1$, preserves positivity and decreases the total number of units. The terminal positive multiple represented with $m$ digits is $p^m-1$ itself. The lower bound on digit sum and (2) give the upper bound on the exponent. At $V=1$ the outer binomial is the weight, a unit, while the inner binomial has valuation $m$ and the predecessor has valuation $k+h$.

On $(1,2^k,2^k)$ one can also count every positive loss threshold. For $1\le j\le k$, put $d=k-j$. When $d=0$ only $V=1$ attains the threshold. When $d>0$ the number of such $V\le X$ is
$$
             \frac{(\log_2 X)^d}{d!}+O_k((1+\log X)^{d-1}).
$$
Indeed (2) requires an odd $V$ with at most $d+1$ ones. Conversely take the units bit and exactly $d$ further ones, all at positions greater than $k+1$ and with pairwise gaps greater than $k+1$. Then $v_2((1+2^{k+1})V-1)=k+1$, the equal-part binomial costs $d+1$ carries, and the outer addition of $V$ and $2^{k+1}V$ has none. This yields defect exactly $j$. Omitting a bounded low position or a bounded-gap pair costs $O_k((1+\log X)^{d-1})$ choices. Rounding the maximum bit length has the same error.

The [checker](check_sparse_predecessor_scales.py) evaluates (1) by factorial floors at 101,979 scales, including every scale through 100,000 and members of both sequences with exponents beyond 1,000. It also checks 12,000 digit-sum bounds and compares the general ray algorithm with the exact exponents above.

~~~sh
python3 check_sparse_predecessor_scales.py
~~~
