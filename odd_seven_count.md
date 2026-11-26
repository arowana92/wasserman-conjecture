# A quadratic count for seven digit units at odd predecessor primes

Among all $N\le X$ for which **some odd prime $p\mid N-1$ has $s_p(N)=7$**, at most $O((1+\log X)^2)$ can fail Wasserman's three-row conjecture. This bound counts weights, including all outer shifts, and is uniform over the choice of the prime. Among normalized exponent patterns of span at most $H$ at these primes, at most $O((1+H)^2)$ can fail at any shift. The earlier [general sparse-position estimate](effective_sparse_positions.md) gives a cubic count at seven units.

At $p=2$ the same quadratic estimates hold for putative exceptions admitting a row avoiding 2 with the units-place 1 in a pair. The other binary allocation has [three gaps $n,n,n\pm1$](seven_binary_alignment.md).

The [large-span theorems](seven_small_prime_span.md) already leave only finitely many possible patterns at all $p\ge5$ combined, with $p\ge29$ excluded altogether. The work here concerns $p=3$ and the stated binary allocation. We give explicit necessary relations for both.

## A half-span bound on the first common divisor

Let $d\ge1$, $0\le e\le f$, $H\ge\max(d,f,8)$, and $p$ prime. If $p=2$, assume $0,e,f$ are distinct. Then
$$ \log_p\gcd(1+p^d,1+p^e+p^f)\le (H+1)/2+6H^{2/3}+\log_p4. \tag{1}$$
The leading coefficient $1/2$ is sharp, for example with $(d,e,f)=(3t,2t,4t)$ at any prime: $p^{2t}-p^t+1$ divides both numbers.

Use the [binomial–trinomial line bound](binomial_trinomial_gcd.md#bounding-the-remaining-scalar). There is one primitive nonnegative direction $(n,u,v)$ with a scalar loss at most $p^{6H^{2/3}}$. Its primitive polynomial gcd has degree $m\le2$. The small internal heights below eight satisfy (1) directly. When the gcd is nonconstant $n>0$, its roots have modulus $p^{-b}$ with $b=d/n$, and $e-bu$, $f-bv$ have absolute value at most one by the same classification. Thus $b\max(n,u,v)\le H+1$. Its constant coefficient is $\pm1$ because it divides the primitive binomial with constant coefficient one. Its value at 1 is bounded by $(1+p^b)^m$.

In fact $2m\le\max(n,u,v)$. Write a common root as $\xi p^{-b}$. If all three heights coincide, primitivity gives $\xi$ of order 6 and the residues $u,v=2,4\bmod6$ in some order. Otherwise the least two heights agree, and the relative magnitudes are $1,1,\sqrt p$ for $p=2,3$, or $1,1,2$ at $p=2$. In the square-root cases, normalize one lowest phase to one. The two lowest phases differ by angle $\pm\pi/2$ or $\pm\pi/3$ respectively, and the opposite phase of their sum has order 8 or 12. Since $\xi^n=-1$, this forces $n\ge4$ or $n\ge6$. In the last binary case the two lowest phases are equal, and the third is opposite. Primitivity and $\xi^n=-1$ force $\xi=-1$, hence $m\le1$. The lowest two slopes differ by a nonzero even integer: equal height and equal slope would give repeated actual positions. This proves the degree claim and (1).

## A second restriction on the six positions

Let $H\ge10^{24}$ be the normalized span of a carry-free pattern with six positive units at $p=2,3$, in the allocation where 1 belongs to a pair. Write the inner blocks as
$$ x=X^i(1+X^d),\qquad y=X^j(1+X^e+X^f),\qquad z=X^k.$$
Let $K=1200000$ at 2 and $K=600000$ at 3. The [common gap scale](sharp_trinomial_gcd.md#consequence-for-seven-unit-counterexamples) gives integers $L\ge H/K$, $|s|,|s'|\le4$ such that
$$ 2d=Ln,\qquad 2e=Lu+s,\qquad 2f=Lv+s'. \tag{2}$$
All three integer slopes are nonnegative, at most $D_0=2K+1$, and $n>0$. In addition, for **at least one** of the translations $t=k-j$ and $t=i-k$, there exist integers $q,r$ such that
$$ 1\le q\le32,\qquad |r|\le D_0 q,\qquad |2qt-rL|<3000000000. \tag{3}$$
The constants are absolute. In particular the residual on the left is a bounded integer, so (3) produces finitely many exact affine equations.

To prove this, the necessary [allocation inequality](seven_span_and_gaps.md#which-allocations-need-checking) reads $27p^{2H}\le U^3A^2C^2$. By (1), $\log_p U\le(H+1)/2+6H^{2/3}+2$. Consequently $A\ge p^{H/9}$ or $C\ge p^{H/9}$: otherwise the exponent on the right would be at most $(3/2+4/9)H+18H^{2/3}+8<2H$.

Put $b=L/2$ and set $F(S)=1+S^n$, $Y(S)=1+p^{s/2}S^u+p^{s'/2}S^v$. At $(S,R)=(p^b,p^t)$, $A$ divides $F(S),Y(S)+R$ for $t=k-j$, while $C$ divides $Y(S),1+R F(S)$ for $t=i-k$. These divisibilities hold away from $p$, including negative $t$. Use whichever divisor is at least $p^{H/9}$ and choose $1\le q\le32$ with $|qt/b-r|\le1/32$ by Dirichlet approximation. Then $|r|\le2Kq+1\le D_0q$.

Substitute $S=T^q$, $R=p^\eta T^r$ with $\eta=t-br/q$. Strip Laurent monomials and clear the least coefficient powers of $p$ to obtain polynomials with algebraic-integer coefficients. At $T=p^{b/q}$ the evaluated pair still has the assigned divisor in the ring with $p$ inverted. The first degree $m_1$ is at most $D_0q$ and $m_1 b\le q(H+2)$. The second degree $m_2$ is at most $2D_0q$. Their coefficient exponent spreads before collecting terms are at most 4 and $4+|\eta|$ respectively. For a nonzero resultant, the absolute value at every embedding is bounded by
$$ 4^{m_1+m_2}p^{4m_2+(4+|\eta|)m_1}. $$
This follows from the Sylvester column sums. Taking the norm bounds the assigned divisor by the same quantity, since it is prime to $p$. Its base-$p$ logarithm would thus be at most $50(K+1)q+(H+2)/32<H/9$. Constant inputs are covered by this argument if the other degree is positive. If both are constant, the first is a normalized $Y$ with at most three terms of exponent spread at most four, so its norm bound suffices directly. Hence the polynomials have a common nonzero root.

At that root write $\alpha=T^q$. In the first case $F(\alpha)=0$ and
$$ p^{qt-br}=(-Y(\alpha))^q\alpha^{-r}; $$
in the second $Y(\alpha)=0$ and
$$ p^{qt-br}=(-1/F(\alpha))^q\alpha^{-r}. $$
The factors displayed in these identities are nonzero. For the first case $\alpha$ is a root of unity, and the absolute logarithmic height of $Y(\alpha)$ is at most $\log(3p^4)$. Indeed its conjugates have modulus at most $3p^2$ and $p^2Y(\alpha)$ is integral.

For the second case $h(\alpha)\le\log(6p^4)\le7\log p$, where $h$ denotes the absolute logarithmic height. To see this, multiply $Y$ by $p^2$ and collect its terms. All coefficients are integral and each archimedean magnitude is at most $3p^4$. At any place, a root satisfies $\max(1,|\alpha|)\le c\max_i|a_i|/|a_{\rm lead}|$ for the collected coefficients, using $c=2$ at archimedean places by the Cauchy bound, $c=1$ elsewhere. The product formula cancels the denominator in the sum of the local logarithms. Also $h(1+\alpha^n)\le n h(\alpha)+\log 2$ by the local triangle inequalities. Using $h(x^{-1})=h(x)$, $h(xy)\le h(x)+h(y)$, and $h(p^z)=|z|\log p$ for rational $z$ now gives $|qt-br|\le6q$ in the first case and $|qt-br|\le7(nq+|r|)+q$ in the second. Doubling gives (3) at both primes. [Silverman's notes, section 4 on heights and the product formula](https://swc-math.github.io/aws/2024/2024SilvermanNotes.pdf)

## Counting patterns and weights

Fix the slopes, offsets and choice of translation in (2) and (3), including the integer value of $2qt-rL$. There are only finitely many such choices independently of $H$. For spans at most $B$, the scale $L$ is a positive integer at most $2d\le2B$. Once it is chosen, the two blocks' internal gaps are fixed, as is one of $k-j,i-k$. The other translation has at most $2B+1$ choices. The condition that the minimum position is zero fixes the common translation. This gives $O((1+B)^2)$ patterns. Bounded spans and the finitely many larger-prime patterns can be included in the constant.

Finally, a fixed normalized pattern of positive span $H$ can fail at at most $315$ shifts in any interval of length $H/9$, by the [explicit shift-spacing theorem](sparse_prime_counts.md#an-explicit-shift-bound-for-seven-units). Put $B=\log_2 X$. For spans in $[R,2R)$ there are $O(R^2)$ patterns and each admits $O(1+B/R)$ shifts up to $B$. Sum over dyadic $R\le B$ to get $O((1+B)^2)$ weights. There are at most finitely many eligible primes, and bounded spans contribute $O(1+B)$. This proves the claimed counts.

Run the following with python-flint installed to check the degree ratios, constants and examples of both exponent-line systems.

~~~sh
python3 check_odd_seven_count.py
~~~
