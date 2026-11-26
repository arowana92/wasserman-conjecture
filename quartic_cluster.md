# Four nearby prime powers

Let $0<a<b<c$ and $P$ be positive integers such that $P,P+a,P+b,P+c$ are powers of four distinct primes. Put $F(X)=X(X+a)(X+b)(X+c)$, $n=F(P)$, and $K=ab+ac+bc$.

* **At $n+1$** the full trinomial conjecture holds if $P>4c^4+1$ and
$$P+b>(a-1)+(b-1)b(b-a)^2(c-b). \tag{1}$$
In particular it suffices that $P>4c^4+1$ and $3125P\ge256c^5$.
* **At both $n$ and $n+1$** it suffices that $P\ge(8c)^{11}$.

Thus a sufficiently small constant times $P^{1/5}$ is an admissible spread at the successor, and a sufficiently small constant times $P^{1/11}$ works at both weights. The interval argument in the [arbitrary-degree theorem](cluster_factor_bounds.md) yields $\gg X^{2/5}/(\log X)^4$ distinct such successors through $X$, and $\gg X^{7/22}/(\log X)^4$ distinct consecutive pairs in the second range.

Under additional conditions on the two inner gaps, the [three-prime certificate](quartic_prime_certificates.md) allows a cube-root spread simultaneously at both weights and also allows a growing multiplier.

## The successor forces a single assignment

Use the [endpoint assignment](endpoint_reduction.md) given by the [all-weight theorem](all_weights_endpoint_reduction.md). Name the row whose coefficient avoids the base of $P$ as $A$ and the companion as $B$. Assign each prime dividing $N(N-1)$ to a row avoiding it. The assigned weight scales are $U,V$ with $UV=N$, and the assigned predecessor factors are $S,T$ with $ST=N-1$. The rows consist of multiples of $U,V$, respectively.

Take $N=n+1$. The base-$P$ blocks of $A$ are three nonnegative coefficient splits of $F(X)+1$ evaluated at $P$. At an assigned node $s\in\{0,a,b,c\}$ their evaluations at $-s$ give the actual unit vector of residues modulo $P+s$ as an equality of integers. In fact the sum of the even-degree nonconstant terms at $-s$ is $s^4+Ks^2$, which equals the absolute sum of the odd-degree terms there. Thus each evaluation lies between $-4c^4$ and $4c^4+1$, and at $s>0$ the strict bound on $P$ rules out any other integer congruent to its required residue; at zero the assertion is immediate. Two such vectors are equal if their arguments differ by at least 2.

At least one side is assigned at least three of the four predecessor powers. Indeed either $S\ge V$ or $T\ge U$, since $(V-1)(U-1)<N-1$. On the side where this inequality holds, each of the three positions must be exceptional at least once. A position exceptional nowhere would, after division by the coprime weight scale, be a positive multiple of the assigned predecessor product and smaller than its normalized row sum. Hence two powers on each side is impossible. If $A$ has exactly three powers it must be this side and the three vectors must be distinct, whereas two of the three nodes differ by at least 2. If $A$ has all four powers then all four vectors coincide: the graph joining integer nodes at distance at least 2 is connected on any four distinct nodes. This again misses a position.

Consequently $A$ receives just $P$. The companion receives $w_1=P+a,\ w_2=P+b,\ w_3=P+c$, with distinct exceptional positions. Write its parts in that order as $V k_i w_jw_l$ for $\{i,j,l\}=\{1,2,3\}$. All $k_i$ are positive integers and $U=\sum_i k_i w_jw_l$.

One part of the binary- or higher-base carry-free row $A$ contains the leading block $P^4$, so the two others sum to at most $N-P^4$. One of these has residue zero at $P$ and equals $UPf$ for some $f\ge1$. It follows that
$$ f\sum_i k_i\le \frac{N-P^4}{P(P+a)(P+b)}<a+b+c.$$
For the strict inequality multiply out $P^3-(P-a-b)(P+a)(P+b)>1/P$. At the block of size $P$ just above the units block, this same part contributes the digit $[Uf]_P$, where brackets denote the least nonnegative residue. The total budget is $abc$. Since $f\sum k_i<3c$ and $3c^3<P$, the residue equals the integer $f(k_1bc+k_2ac+k_3ab)$. Thus $f\sum_i k_i/h_i\le1$ for $(h_1,h_2,h_3)=(a,b,c)$. In particular $k_1\le a-1$, $k_2\le b-1$.

## The carry at the middle power

At $Y=w_2$ the part $B_1$ is divisible by $Y$. Its next digit $r=[B_1/Y]_Y$ is at most the corresponding digit of $N$, by carry freedom at the base prime of $Y$. That digit is
$$[P w_1 w_3]_Y=b(b-a)(c-b),$$
using $P>4c^4+1$. Also $B_2=Vk_2w_1w_3\equiv1\bmod Y$ and $r\equiv V k_1 w_3\bmod Y$. Multiplying the latter congruence by $k_2w_1$ gives
$$0<k_1+k_2(b-a)r\equiv0\pmod Y.$$
This forces $P+b\le(a-1)+(b-1)b(b-a)^2(c-b)$, proving (1). For the uniform estimate use $b^2(b-a)^2(c-b)\le b^4(c-b)\le256c^5/3125$.

The same digit argument works at an arbitrary predecessor $m w_1w_2w_3$ when the three powers are assigned to distinct exceptional positions in one row. For $w_1<Y=w_2$, the general inequality is $Y\le k_1+k_2(Y-w_1)[m w_1w_3]_Y$. It uses no spacing assumption.

## The product

Take $N=n$, $M=n-1$, and suppose $P\ge(8c)^{11}$. As in the [general proof](cluster_factor_bounds.md), the parts of $A$ are polynomials $B_i(P)$ splitting the coefficients of $F$. For every assigned weight power $P+s$ all part polynomials vanish at $-s$ by the evaluation bound above. Their product of assigned linear factors $\Phi$ divides each part polynomial. At most two powers are assigned to $A$: the sign count of the general proof gives $m\le2$ at degree four. Its digit hypotheses hold here, since every translated coefficient has absolute value at most $(2c+1)^4<P/4-1$. Write $L_i=B_i/\Phi$, so $L_i(P)$ are the parts normalized by $U=\Phi(P)$. Each full power of $S$ divides at least two of them.

The polynomial $F-1$ is irreducible, with roots $|\alpha|\le c+1$, as proved in the general argument. If $m=1$ then $\Phi=X$ and each $|L_i(\alpha)|\le(2c+1)^3$. Thus $S^2\le(2c+1)^{36}$ by taking the three resultants with $F-1$. But $27T^2\le P^3$, so $S>P^{5/2}$ and $P^5<(2c+1)^{36}$, contrary to the threshold.

Suppose $m=2$ and $\Phi=X(X+s)$. The other two weight powers have product $V=(P+t)(P+u)$. Two part polynomials have degree at most three since their leading coefficients sum to one. After division by $\Phi$ these two are $x=\alpha X+\beta$ and $y=\gamma X+\eta$, the third being $V(X)-x-y$, with $V(X)=(X+t)(X+u)$. The linear slopes lie in $[0,3c]$, and $|\beta|,|\eta|\le3c^2$, by comparing the cubic and quadratic coefficients before division. The resultant bound $|\operatorname{Res}(F-1,x)|\le(7c^2)^4$ and its counterpart for $y$ follow from the root bound and $c\ge3$. Every power of $S$ divides $x(P)$ or $y(P)$.

We claim $S<(8c)^{11}$. If one polynomial is a constant, that constant is positive and at most $3c^2$ and the resultant bound for the other suffices. If the two nonconstant polynomials are proportional, write them as positive integer multiples of one primitive linear polynomial. The multipliers are at most $3c$, and $S$ divides their product times the gcd of $M$ with the primitive polynomial value. The latter is at most $(7c^2)^4$ by the same root bound.

In the remaining case, powers dividing both $x(P),y(P)$ have product at most $|\alpha\eta-\beta\gamma|\le18c^3$. The others divide either the pair $x(P),V(P)-y(P)$ or the pair $y(P),V(P)-x(P)$. The corresponding linear–quadratic resultants are at most $72c^4$ in absolute value. For example the first is
$$\beta^2-(t+u-\gamma)\alpha\beta+(tu-\eta)\alpha^2.$$
If neither resultant vanishes, $18c^3(72c^4)^2$ bounds $S$.

If one vanishes, say the first, the root $z=-\beta/\alpha$ is an integer root of the monic quadratic $V-y$. Its coefficients give $|z|^2\le3c|z|+4c^2$, so $|z|\le4c$. All powers dividing $x(P)$ together contribute at most $\alpha|F(z)-1|\le1878c^5$, since $P-z$ divides $F(P)-F(z)$ and $F(z)-1\ne0$. The powers left divide $y(P),V(P)-x(P)$. Their product is at most $72c^4$ if the other resultant is nonzero, and at most $1878c^5$ by the same root argument otherwise. These bounds all give the claim. Finally $U=P(P+s)<2P^2$ and $27T^2\le U^3$, so $S=M/T>P$. This contradicts the claim and proves the product assertion.

## An explicit pair

Let $P=10^{45}+9$. The four numbers $P,P+112,P+292,P+312$ are prime. Their product and its successor both have 181 decimal digits and satisfy the conjecture. Run `python3 check_quartic_cluster.py` with python-flint installed to verify the primalities by recursive order certificates, the numerical thresholds, and the carry and polynomial calculations used here.
