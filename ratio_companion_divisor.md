# A ratio in one row and an arbitrary companion

Let $A=(ag,bg,z)$ and $B$ be positive trinomial rows of weight $N$, with $\gcd(a,b)=1$. Write $B$ as an integer multiple of its primitive direction $y=(y_1,y_2,y_3)$, and put
$$
 u=y_1+y_2+y_3,\qquad d_i=\gcd(y_j,y_k)\quad(\{i,j,k\}=\{1,2,3\}),\qquad H=d_1d_2d_3 .
$$
The $d_i$ are pairwise coprime. Let $C=C(a,b)$ be the [sharp ratio constant](part_ratio_predecessor.md), $M=N-1$, and $G=\gcd(NM,K_A,K_B)$. Then
$$
 uM\mid Gg\,\gcd(M,\operatorname{lcm}(C,H)),\qquad
 G\ \ge\ \frac{u(a+b)}{\gcd(M,\operatorname{lcm}(C,H))}.               \tag{1}
$$
If $G=1$ the denominator in the second formula strictly exceeds $u(a+b)$. In particular, **one integer or geometric ratio suffices against every primitive companion with $H\le u$**, without any constraint on the scale or the third part in $A$.

To prove (1), first $u\mid K_B$ because $N$ divides each part of $B$ times $K_B$, and $y$ is primitive. Also $N\mid gK_A$ by the same factorial identity. Therefore $u\mid g\gcd(N,K_A,K_B)$. At each $p^e\parallel M$ let $w=p^{e-\min(e,v_p(K_A),v_p(K_B))}$. When $w>1$, the row of smaller valuation has residues $0,0,1$ modulo $w$ by the [retained-power lemma](endpoint_pair_bound.md). If this is $B$, then $w\mid H$ since its scale divides $N$. If this is $A$ and $p\mid g$, the first two parts both vanish, so $w\mid g$. If this is $A$ and $p\nmid g$, the defining divisibility $M\mid gCK_A$ gives $w\mid C$. Thus $M/\gcd(M,K_A,K_B)$ divides $\operatorname{lcm}(g,C,H)$. Combining the two statements prime by prime proves the divisibility in (1). The inequality uses $(a+b)g\le M$. When $G=1$ one has $u\mid g,N$ and hence $u\mid z$. As $u\ge3$, the last inequality is strict.

## A full homogeneous relation

Suppose $r y_1+s y_2=t y_3$ with $r,s,t$ positive and $\gcd(r,s,t)=1$. Then
$$
 \frac{u}{H}\ \ge\ \frac1r+\frac1s+\frac2t.                            \tag{2}
$$
Indeed $d_1\mid r$, $d_2\mid s$ and $d_3\mid t$, by primitivity. Each $y_i$ is at least $H/d_i$. Eliminate $y_3$ to get $u=(1+r/t)y_1+(1+s/t)y_2$. Substituting the lower bounds gives (2). It is sharp on every pairwise coprime set of coefficients $r,s,t$. Take the primitive direction proportional to $(st,rt,2rs)$. Before normalization its pair gcds are $r\gcd(t,2)$, $s\gcd(t,2)$ and $t$. Its coordinate gcd is $\gcd(t,2)$. Substitution gives equality.

Consequently **the common prime exists at every weight** whenever
$$
                \frac{C(a,b)}{a+b}\le\frac1r+\frac1s+\frac2t.           \tag{3}
$$
A coefficient of magnitude one anywhere in the full relation suffices against every ratio $1:b$ and every ratio $(1+P+\cdots+P^{m-1}):P^m$ with $P$ a prime power. If $t=1$, every reduced ratio with smaller term at most two works, since $C(a,b)\le ab$. Against $y_1+y_2=y_3$, every reduced ratio with smaller term at most four works. More generally (1) gives $G\ge(a+b)(1/r+1/s+2/t)/C(a,b)$.

If two rows each satisfy a full relation with a coefficient of magnitude one, they also share an endpoint prime at every weight. To see this, assume $G=1$. Their primitive sums $u,v$ are coprime since each divides its own coefficient. Hence $uv\mid N$. The retained powers give $M\mid\operatorname{lcm}(H_y,H_{y'})$, whereas (2) gives $H_y<u$ and $H_{y'}<v$. Their product is less than $uv-1\le M$, a contradiction.

## How many companion directions escape the bound?

Fix $Y\ge1$. Among all ordered primitive triples of sum at most $X$, all but $O(Y X^2\log^2(2X))$ work **simultaneously at every scale against every ratio** with $C(a,b)/(a+b)\le Y$. This includes every reduced ratio with smaller term at most $Y$. For any fixed $Y>0$, the number of directions with $H>u/Y$ has exact order $\Theta_Y(X^2\log^2 X)$. This last count concerns the numerical test on a single direction. Actual counterexamples would have to satisfy the full carry conditions as well.

For the upper bound restrict $R/2<u\le R$. At given $H$ there are at most $\tau_3(H)$ ordered factorizations $H=d_1d_2d_3$. For each, $H/d_i\mid y_i$ leaves at most $R^3/H^2$ triples. Since $\sum_{n\le L}\tau_3(n)\ll L\log^2(2L)$ by summing over the first two factors, dyadic summation gives $\sum_{n>T}\tau_3(n)/n^2\ll\log^2(2T)/T$ for $T\ge1$. Use $T=R/(2Y)$, or the trivial $O(R^3)$ when $T<1$, and sum over $R$. Equation (1) handles every direction outside this tail.

For the matching lower bound choose dyadic starting values $D_1,D_2$ independently in $[X^{0.30},X^{0.32}]$. Take $D_3$ a power of two in $[X/(D_1D_2),2X/(D_1D_2))$. Take six pairwise coprime integers
$$
 D_i\le d_i<2D_i,\qquad \theta D_i\le z_i<2\theta D_i,\qquad
 \theta=\min(1,Y)/120 .
$$
Each box supplies $\gg_Y(D_1D_2D_3)^2\gg_Y X^2$ choices. Indeed, excluding joint divisibility of two coordinates by each prime up to a fixed $P$ leaves an asymptotic proportion $\prod_{p\le P}((1-1/p)^6+(6/p)(1-1/p)^5)$ by the CRT. This product has a positive limit. Primes above $P$ spoil at most $O_Y(\sum_{n>P}n^{-2})$ of the full box, since the count of multiples of $p$ in each interval is bounded by its upper endpoint divided by $p$. Choose $P$ sufficiently large, then let $X$ grow. The least side length tends to infinity uniformly.

Put $H=d_1d_2d_3$ and $y_i=z_iH/d_i$. The pair gcds recover exactly the $d_i$ and the triple is primitive. Its sum is at most $6\theta H$, with $H<16X$, so $u\le X$ and $H>u/Y$. Different choices of the six integers give different ordered directions, and there are $\gg(\log X)^2$ disjoint choices of the first two dyadic intervals.

The [checker](check_ratio_companion_bound.py) verifies the divisibility and inequality in (1) on 856,348 exact coefficient comparisons, checks (2) on 60,856 plane incidences and verifies 580 sharp examples. The numbers of unordered primitive directions with $H>u$ through sums 32, 64, 128 and 256 are respectively 0, 6, 63 and 508.

~~~sh
python3 check_ratio_companion_bound.py
~~~
