# Seven digit units at predecessor primes 5 and 7

Let $p\in\{5,7\}$, $a\ge1$ and $0=h_1\le\cdots\le h_6=H$. For every choice of the six positions and every outer shift $a$, Wasserman's three-row conjecture holds at
$$N=1+\sum_{j=1}^6 p^{a+h_j}\qquad\text{if }H\ge10^{16}.$$

Together with the [large-span theorem at $p\ge11$](seven_large_span.md) and the [all-span theorem at $p\ge29$](seven_base27.md), this leaves only finitely many normalized seven-unit patterns at **all predecessor primes $p\ge5$ combined** that could fail at any shift. The bound is on the difference between the highest and lowest positive exponents. It places no upper bound on the first exponent or on the weight.

## The three common factors

If any specified terms carry, the base-$p$ digit sum of $N$ drops to at most five and the [five-unit theorem](six_digit_gaps.md) applies. In the remaining case the [Kummer reduction for seven units](seven_span_and_gaps.md#which-allocations-need-checking) gives the following necessary conditions for a counterexample. Put
$$x=X^i F,\quad F=1+X^d,\quad y=X^j(1+X^e+X^f),\quad z=X^s,$$
where the six positions are those of $x,y,z$, $d>0$ and $0\le e\le f$, $f>0$. There must exist odd integers $U,A,C$, prime to $p$, such that
$$U\mid F(p),y(p),\qquad A\mid F(p),(y+z)(p),\qquad C\mid y(p),(x+z)(p),\qquad 27p^{2H}\le U^3 A^2 C^2. \tag{1}$$
Here $U\ge7$. Moreover $\gcd(C,p^{3d}-1)=1$. Indeed $U$ divides both $p^d+1$ and $p^{a+s}+1$, so $v_2(d)=v_2(a+s)$ by taking the order modulo any prime divisor of $U$. Each prime of $C$ divides $p^{a+s}+1$ too, hence its order for $p$ cannot divide $3d$. The [trinomial factor argument](seven_large_bases.md#bounding-the-factor-from-two-trinomials) shows that after stripping monomials the polynomial gcd of $y,x+z$ divides $X^{3d}-1$. In particular its value is coprime to $C$.

Set $h_U=\gcd(F,y)$, $h_A=\gcd(F,y+z)$ monic, and $k=\deg h_A$. Four unit complex numbers summing to zero have an opposite pairing. Each matching contributes, if nonempty, the roots of $X^{d/n}+1$ for an odd divisor $n\ge3$ of $d$. There are at most two such matchings. For two distinct indices $n_1,n_2$ the fraction $k/d$ equals $1/n_1+1/n_2-1/\operatorname{lcm}(n_1,n_2)$. In particular $k\le7d/15$. The case $F\mid y+z$ is impossible since it would give $U\mid N-1,N$. The monic lcm of the matching factors gives $h_A(p)\le4p^k$ (divide their product by the common binomial if there are two). See the [matching argument](seven_large_bases.md#common-roots-at-the-unshifted-pair).

Apply the [sparse binomial gcd estimate](binomial_sparse_gcd.md) after stripping the lowest monomial. It gives
$$U\le h_U(p)p^{6H^{2/3}},\qquad A\le 4p^{k+7H^{3/4}}.$$
Use the small-height trivial bound if a block has internal height below the threshold of the estimate. In the bounds below the combined logarithmic error from $U^3A^2$ is at most $E=18H^{2/3}+14H^{3/4}+2\log_p4$. For $H\ge10^{16}$,
$$ E+2\log_p8+2 < H/500. \tag{2}$$
Indeed $H^{1/3}\ge200000$, $H^{1/4}\ge10000$ and $2\log_5 4+2\log_5 8+2<8$. All logarithmic comparisons below can be checked with integer powers.

## When the binomial and trinomial have a common root

If $h_U\ne1$, write $t=\gcd(d,e,f)$. The [three-term root classification](binomial_trinomial_gcd.md#locating-every-common-root) gives $h_U=X^{2t}-X^t+1$, $d/t\equiv3\pmod6$, and $e/t,f/t$ are $2,4\bmod6$ in some order. Thus $H\ge4t$ and no difference in $y$ has the 2-adic valuation of $d$. An opposite pairing of $y+z$ would pair two terms of $y$, hence $k=0$. Since $y$ has distinct positions the [companion bound](seven_base31.md#two-bounds-on-the-companion-factor) gives $C\le3^H$. If $H\ge10t$, then (1) fails by $6t/H+2\log_p3+1/500\le3/5+2\log_53+1/500<2$.

For $H<10t$ the two positive gap quotients in $y$ are $(2,4)$ or $(4,8)$. In either case $1+X^e+X^f$ divides $X^{12t}-1$. Apply the sparse binomial estimate with the minus sign to $X^{12t}-1$ and $x+z$ stripped of its lowest monomial. Their polynomial gcd divides $X^{3d}-1$, since any root of unity at which the three terms of $x+z$ sum to zero has its $d$-th power a primitive cube root and the binomial is squarefree. Its value is prime to $C$. Thus $\log_p C\le6(3H)^{2/3}<13 H^{2/3}$, as $12t\le3H$. Now $6t\le3H/2$ and $3H/2+E+26H^{2/3}<2H$, again contradicting (1).

## When the four-term polynomial must supply a matching

Assume henceforth $h_U=1$. If $C>3^H$, the [repeated-position bound](seven_base31.md#two-bounds-on-the-companion-factor) gives $k\le1$ and $C\le4^H$. This is impossible since $2\log_54+1/500<2$. So $C\le3^H$. Equations (1) and (2) now require
$$k/H>11/35,\qquad k/H>3/7\quad\text{if }p=7. \tag{3}$$
Indeed $22/35+2\log_53+1/500<2$ and $6/7+2\log_73+1/500<2$. Without index 3 the largest possible matching fraction $k/d$ is $11/35$, attained by 5 and 7. The other pairs give at most $1/5+1/9$. Consequently $d=3w$ and the four terms of $y+z$ pair with gaps odd multiples of $w$. Since $k\le7w/5$, (3) gives $H<49w/11<5w$. Both gaps therefore equal $w$ or $3w$.

### Repeated positions in one of the two trinomials

We claim that if either $x+z$ or $y$ repeats a position, then $C\le 8\cdot2^{H+w/3}\le8\cdot2^{10H/9}$. If $x+z$ repeats, its gap is $3w$ and $y$ contains a complete matched pair at distance $w$ or $3w$. If $y$ repeats, its two equal terms must both be opposite to the other two at a primitive $2w$-th root in the index-3 matching. All four positions of $y+z$ are congruent modulo $w$. Also the two matchings obtained by interchanging the equal terms have identical root sets, so $k=w$ and $H<35w/11$. The nonzero gap of $y$ is $nw$ with $1\le n\le3$, while $x+z$ contains the gap $3w$.

Strip monomials and reverse both inputs if required to put the repeated-position polynomial in the form $X^{nw}+2$. Reversal preserves absolute resultants. Relative to the start of the pair in the other input, its support is $(0,Lw,c)$ with signed $c$. Here $(n,L)$ is $(3,1),(3,3),(1,3)$ or $(2,3)$. The full span gives $\max(c,Lw)-\min(c,0)\le H$ and $3w\le H$. Put $g=\gcd(w,c)$, $m=w/g$, $b=c/g$. Before reversal the integer $C$ divides the resultant of the two polynomials written in the common variable $X^g$, by evaluation at $p^g$. The resultant is nonzero: the repeated binomial is Eisenstein or its reciprocal and the other polynomial either has only algebraic-unit roots or is another such binomial. Identical stripped forms are impossible since $U$ would then divide $z(p)$.

For a root $S$ of $S^n=-2$, put $\rho=|S|=2^{1/n}$. Over the $m$ roots of $T^m=S$ the contribution to the absolute resultant is bounded by
$$\begin{cases}
 |1+S^L|^m+\rho^b & b\ge0,\\
 1+\rho^{-b}|1+S^L|^m & b<0.
\end{cases}$$
Indeed $\gcd(b,m)=1$ permutes the $m$ phases in $T^b$, and for $b<0$ the stripped polynomial has the additional monomial $T^{-b}$. Each sum is at most twice its larger term.

For $L=3$, the magnitudes $|1+S^3|$ are respectively $7,3,1$ at $n=1,2,3$. Each is at most $\rho^{3+1/3}$. For $b\ge0$ use $c,3w\le H$ and for $b<0$ use $|c|+3w\le H$. The larger term is at most $\rho^{(H+w/3)/g}$ in each factor. For $(n,L)=(3,1)$ the two nonreal roots have $|1+S|=\sqrt{1+\rho+\rho^2}<2=\rho^3$. The negative real root gives $|1+S|=\rho-1<\rho^{-2}$. When $b\ge0$ all three maxima are at most $\rho^{H/g}$. When $b<0$ we have $-b\le H/g-m$. The real maximum is at most $\rho^{H/g-3m}$, and the other two are at most $\rho^{H/g+2m}$. Multiplying proves the claimed bound, in fact $C\le8\cdot2^{(H+w/3)/g}$. This contradicts (1) since $14/15+(20/9)\log_52+1/500<2$.

### Distinct positions in both trinomials

In this case $C\le 8^{H/2}$. Use the [Sylvester column minor for the primitive gcd](seven_large_bases.md#bounding-the-factor-from-two-trinomials), whose value at $p$ is prime to $C$. Each trinomial contributes at most $H$ shifts. Pad its collection of possible shifts to the indices $0,\ldots,H-1$. Of the two consecutive support gaps choose the smaller, $l\le H/2$. Columns whose shifts differ by $l$ have squared norms 3 and inner product at least 1. In each residue class modulo $l$ there are at least two consecutive shifts along the chain. Partition every chain into consecutive pairs and triples.

For a full segment of length $r=2,3$, orthogonalization bounds its volume by $\sqrt3(8/3)^{(r-1)/2}\le8^{r/4}$, using the projection on the preceding column. A proper subset has volume at most $3^{(r-1)/2}\le8^{r/4}$. Projection onto minor rows can only decrease volume, and the volume of all columns is at most the product over the segments. Using both shift collections proves the bound. Since $7/10+\log_58+1/500<2$, (1) now requires $k/H>7/20$.

In particular $k>d/3$, so there are two matching indices $3,n$, where $n\ge5$ is odd and coprime to 3. Put $g=d/(3n)$. Then $k=(n+2)g$. The two matchings form a connected graph on the four positions and every matched difference is divisible by $g$, so all four positions agree modulo $g$. From $7d/20<k$ we get $n<40$, and $H/g<20(n+2)/7<120$.

## Bounding the companion factor when the matching degree is large

Write $y=X^j Y(X^g)$, where $Y$ is a monic trinomial with constant 1. Put $v=i-s$. At $X=p$, $C$ divides both $Y(p^g)$ and $1+p^v(1+p^{3ng})$ away from $p$. We first show that $Y(S)$ and $1+S^{3n}$ have resultant prime to $p$.

At $p=5$, the complete pair in $y$ from the index-3 matching has gap $mn$ after division by $g$, for $m=1$ or 3. After Laurent monomial normalization, $Y$ has terms at $0,mn,b$ for a signed $b$. At a common root $\alpha$ in characteristic 5 put $\eta=\alpha^n$. We have $\eta^3=-1$. If $\eta=-1$ or $m=3$, the complete pair cancels and the singleton cannot vanish. Otherwise $\eta$ is a primitive sixth root and $\eta^5=\eta^{-1}$. In $\mathbb F_{25}$ its norm is 1, whereas the norm of $1+\eta$ is $2+\eta+\eta^{-1}=3$. Raise $\alpha^b=-(1+\eta)$ to the $n$-th power before taking norms. This gives $1=3^n$ in $\mathbb F_5$, impossible for odd $n$.

At $p=7$, (3) and $k/d=(n+2)/(3n)$ force $n=5$. Now $H/g<49/3$. Subtract the least exponent of the four and divide by $g$. Among $(0,a,b,c)$ with $1\le a<b<c\le16$, only $(0,3,5,8)$ and $(0,5,9,14)$ give degree seven for the gcd of the four-term sum with $1+S^{15}$. This small enumeration checks all $\binom{16}{3}=560$ sets by polynomial Euclidean gcd. Deleting any singleton in either set and stripping the lowest monomial gives absolute resultant 271, coprime to 7. The verifier checks the eight resultants by a Sylvester determinant and independently by SymPy.

Choose $Q=\lceil\sqrt H\rceil$. Dirichlet approximation gives $1\le q\le Q$ and an integer $r$ with $|qv/g-r|\le1/Q$, hence $|r|\le120q+1$. On the line $S=p^g T^q, R=p^v T^r$ consider $Y(S)$ and $1+R(1+S^{3n})$. A common nonzero root would give $\alpha=S$ an algebraic unit and
$$1=(-1)^q p^{qv-gr}\alpha^r(1+\alpha^{3n})^q.$$
The norm of $1+\alpha^{3n}$ in $\mathbb Q(\alpha)$ is an integer prime to $p$, since it divides the resultant just considered. Taking norms and then $p$-valuations forces $qv=gr$.

If $qv\ne gr$ the resultant is nonzero after removing Laurent monomials. It is divisible by $C$ away from $p$ by evaluation at $T=1$. The degrees are at most $120q$ and $240q+1$. In the resultant every term has the same slope weight, and the sum of absolute multiplicities is at most $3^{360q+1}$ by the column sums. Only monomials involving $R$ have residual exponent $v-gr/q$, of absolute value at most $g/(qQ)$. The two terms containing $R$ have the same residual, and at most $120q$ columns use the second polynomial. Thus the spread in exponents of $p$ is at most $120H/Q$. Clearing the lowest power gives $\log_p C\le360q+1+120H/Q\le850\sqrt H$.

If $qv=gr$, reduce the direction to $(q_0,r_0)$ primitive with $q_0>0$. Then $g=bq_0$, $v=br_0$ for an integer $b>0$. Use the unweighted polynomials $Y(T^{q_0})$ and $1+T^{r_0}(1+T^{3nq_0})$, stripping Laurent monomials. Their primitive gcd at $T=p^b$ is prime to $C$: substituting $T=X^b$ in that gcd gives a primitive divisor of the original stripped pair $y,x+z$. The same Sylvester minor argument with column norms at most three now gives $C\le3^{360q+1}$. This proves the square-root bound also on these lines.

Finally (1) would require $2H\le(14/15+1/500)H+1700\sqrt H$, which fails for $H\ge10^{16}$. This proves the theorem.

## Reproducing the computations

Run the following with python-flint and SymPy installed. The checks include 3,885 repeated-position resultants, 26,565 distinct-position resultants, and 332 modular gcds over all 83 eligible normalized four-position profiles at 5. The only finite computation required in the proof beyond rational inequalities is the eight-resultant check at 7.

~~~sh
python3 check_binomial_sparse_gcd.py
python3 check_seven_small_prime_span.py
~~~
