# Intervals at offsets from powers of a prime

Fix $\eta>0$ and $\alpha<6/7$. For all sufficiently large $N$ write
$$N=Q+h,\qquad h=dP+c,\qquad Q=p^k,\quad P=p^j<Q,$$
where $p$ is prime, $d\ge1$ is an integer and $0\le c<P$. Suppose
$$P\ge N^{1/2+\eta},\quad h\le N^\alpha,\quad c\le Q^{2/5}.$$
Then every counterexample must contain $E=(1,1,N-2)$. If also $p\mid N(N-1)$, the full conjecture holds at $N$ under **either** of the following conditions:
$$h^7\le6^6(N-1)^4,\qquad\text{or}\qquad
 h^4\max(1,c)\le (N-1)^2P.$$
The threshold depends only on $\alpha,\eta$. In base two the condition $p\mid N(N-1)$ always holds.

**Proof.** In a counterexample choose a row $A$ avoiding $p$. For large $N$, $h<Q$. Kummer's criterion assigns the leading $Q$ to a single part of $A$. Write its other two parts as $a=iP+u$, $b=j'P+v$, with $i,j'\ge0$, $0\le u,v<P$. Adding the low digits gives $u+v\le c$. Also $s=a+b\le h$.

If $E$ is absent, the [uniqueness theorem](near_row.md) forces $A$ to be the minimizing row. Indeed, choose $\delta>0$ small enough that $\alpha<6/7-\delta$ and $\delta<\eta$. Then $s\le N^{1/2+\delta}<P$. Thus $i=j'=0$ and $s\le c$, contradicting the [lower bound $s>2.73N^{2/5}$](critical_small_parts.md).

If $E$ occurs and $p\mid N(N-1)$, $A$ differs from $E$. By the [two-part endpoint inequalities](endpoint_pair_bound.md), $s^7>6^6(N-1)^4$, excluding the first condition. For the second, take the assignment of factors of $N,N-1$ from that proof, with $U$ the product from $N$ assigned to $A$. Thus $U\mid a,b$ and
$$27(N-1)^2\le ab(s-1)U,\qquad ab\ge3\sqrt3(N-1)\sqrt{\gcd(a,b)}.$$
If $i=j'=0$ the same contradiction $s\le c$ applies. Otherwise $U\mid ib-j'a=iv-j'u$. When this determinant is nonzero, $U\le hc/P$ since each of its two nonnegative terms is at most $hc/P$. Hence $27(N-1)^2<h^4c/(4P)$, contrary to the condition. When it is zero, divide $(i,j')$ by its gcd $t$ to obtain a primitive integer vector. Then $(u,v)$ is a nonnegative integer multiple $\ell$ of the same vector, so $\gcd(a,b)=tP+\ell\ge P$. Since $ab\le h^2/4$, the second displayed inequality gives $h^4\ge432(N-1)^2P$, also impossible.

## Sizes of the certified sets

For each fixed $\epsilon>0$ and all sufficiently large powers $Q=2^k$, at least $\gg_\epsilon Q^{33/70-\epsilon}$ distinct weights in $[Q,Q+Q^{4/7}]$ satisfy the full conjecture by this theorem. To see this, fix $0<\theta<\min(\epsilon,1/14)$ and take $P$ to be the least power of two at least $Q^{1/2+\theta}$. Take
$$1\le d\le Q^{4/7}/(2P),\qquad 0\le c\le Q^{2/5}.$$
The intervals are disjoint since $c<P$. For all large $Q$ the hypotheses hold with $\eta=\theta/2$ and $\alpha=4/7$, and the first endpoint condition applies. The number of weights is $\gg Q^{4/7-1/2-\theta+2/5}$.

The second endpoint condition also gives intervals beyond distance $Q^{4/7}$. For instance take $P=2^{\lceil3k/5\rceil}$, $1\le d\le Q^{1/100}$ and $0\le c\le Q^{1/10}$. Here $h\ll Q^{61/100}$ and $h^4\max(1,c)\ll Q^{254/100}$, while $(N-1)^2P\ge Q^{260/100}$. Thus every weight in every one of these intervals satisfies the conjecture for all sufficiently large $k$.
