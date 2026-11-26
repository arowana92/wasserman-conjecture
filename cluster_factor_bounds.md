# Growing intervals containing products of nearby prime powers

Let $n\ge2$, $0=h_0<h_1<\cdots<h_{n-1}=h$, and $d\ge1$ be integers. Suppose $P+h_i$ are powers of primes. Write
$$F(X)=\prod_{i=0}^{n-1}(X+h_i),\qquad H=h+1,\qquad C=(4^n d)^{6n+2}.$$
The full trinomial common-divisor conjecture holds at the following weights.

| Weight | Sufficient lower bound on $P$ |
|---|---|
| $dF(P)+1$ | $C H^{3n/2}$ |
| $dF(P)$ | $C H^{3n(n-b)/(2n-3b)}$, where $b=\lfloor n/2\rfloor+\mathbf1_{n\equiv1\bmod4}$ |

The constants allow the multiplier to vary. At fixed $n,d$ the exponent of the spread $H$ for the successor is $3n/2$. For the product it is $3n$ at even $n$. The [earlier translated-product theorem](translated_products.md) also handles bounded cofactors of the translates and gives a larger sufficient threshold when all translates are powers. The [four-factor estimates](quartic_cluster.md) sharpen the two exponents to 5 and 11 at $n=4,d=1$.

## Polynomial factors from carry avoidance

Suppose there is a counterexample at $N=dF(P)+\epsilon$, $\epsilon\in\{0,1\}$. Use the [all-weight endpoint reduction](all_weights_endpoint_reduction.md) and the [full prime-power assignment](endpoint_reduction.md). Denote the two nonendpoint rows by $A,C'$, choosing $A$ to avoid the base of $P$. Write $UV=N$, $ST=N-1$ with primes of $US$ assigned to $A$. Its parts are multiples of $U$, while the parts of $C'$ divided by $V$ are positive integers of sum $U$ and product divisible by $T^2$. In particular $U,V\ge3$ and $27T^2\le U^3$.

At the displayed thresholds $h<P$ and the $n$ translated powers have distinct bases. The coefficients of $dF+\epsilon$ are smaller than $P$. Grouping the base-prime digits into blocks of size $P$, Kummer's criterion expresses the parts of $A$ as $B_i(P)$ for three polynomials with nonnegative integer coefficients summing coefficientwise to $dF+\epsilon$. Each is nonconstant since a positive constant could only be 1.

Let $I$ index the translated powers assigned to $A$, put $m=|I|\ge1$, $k=n-m$, and $\Phi(X)=\prod_{j\in I}(X+h_j)$. At each assigned power the three residues are all zero for $\epsilon=0$, or a coordinate unit vector for $\epsilon=1$. These are the **integer values** at $X=-h_j$ as well. Indeed every coefficient of $B_i(Y-h_j)$ in absolute value is at most $d(2H)^n+1<P/4-1$, by summing the absolute coefficients and using the nonnegative coefficient bound. This also bounds the translated coefficients of the total. At the successor the node vectors agree at any two nodes whose difference is at least 2, by divisibility of polynomial value differences. At most two distinct vectors occur, since three integer nodes cannot have all pairwise differences one. Thus some part vanishes at every assigned node. At the product all three parts do. Each such part has the factorization $B_i=\Phi L_i$ over $\mathbb Z[X]$, with $\deg L_i\le k$.

The exclusion of constant parts uses the fact that each part of this row is a positive multiple of the weight scale $U\ge3$.

## A bound after dividing out the assigned factors

Put $G_\pm=dF\pm1$ and $K=(2\cdot3^n)^n d^{2n}$. Every root of $G_\pm$ has modulus at most $H$, since outside that disk every $|X+h_j|>1$. The polynomial $G_-$ is irreducible and $G_+$ is irreducible or the square of an irreducible integer polynomial. In fact both are primitive. In a nontrivial integer factorization $JK$, the values of the two factors at the $n$ nodes are units, with opposite signs for $G_-$ and equal signs for $G_+$. Thus $J+K$ or $J-K$ vanishes identically by its degree. A negative square has the wrong leading coefficient, and a reducible square root would also yield a factorization with unequal factors.

On $|X|=2H$, $|B_i(X)|\le 2d(3H)^n$ and each factor $|X+h_j|\ge H$. Hence $|L_i(X)|\le 2d\,3^n H^k$ throughout the disk by the maximum modulus principle. For an irreducible $G_\pm$, the nonzero integer resultant with $L_i$ has absolute value at most $K H^{nk}$, by multiplying over its roots and including the leading-coefficient factor. The gcd of the two evaluations at $P$ divides this resultant by its integral Bézout identity.

At the successor this bounds $U$, since $U\mid N$ and $U\mid L_i(P)$ by coprimality to $\Phi(P)$. If $G_+=J^2$, write $L_i=J^a L'$ with $a=0$ or $1$, as $\deg L_i<n$. On the same circle $|J(X)|\ge H^{n/2}$, by its leading coefficient and the bound on its roots. Thus $|\operatorname{Res}(J,L')|^2\le K H^{nk}$ by the same argument. For $p^v\parallel J(P)$ and $p\mid U$, unitarity gives $p^{2v}\mid U$, and $p^v\mid L'(P)$. Therefore $U$ divides this squared resultant. This proves the same bound in the square case.

## Counting the missing powers

At the successor $k\ne0$. Otherwise all nodes are assigned. At the largest translated power $w=P+h$ use the variable $Y=X+h$. At the $n-1$ earlier nodes, now positive real arguments, two part polynomials vanish at each. A polynomial with positive leading coefficient and $r$ distinct positive roots has at least $\lceil r/2\rceil$ strictly negative coefficients, by Descartes' rule of signs. The three parts therefore need at least $n-1$ such coefficients altogether. For $n=2$ they need at least two, since two separate parts vanish at the earlier node.

In base $w$, a strictly negative translated coefficient below the leading term gives a digit greater than $w/2$. Normalization of the signed coefficients borrows at most one at a time, by their absolute bounds. Two parts cannot use the same such position. Nor can it occur opposite a strictly positive coefficient of the total, which gives a digit smaller than $w/2$. The constants at the largest node are nonnegative. The other coefficients of $dF(Y-h)+1$ alternate strictly in sign, giving just $\lfloor n/2\rfloor$ available negative positions. This is a contradiction.

Now $P^k\le T\le U^{3/2}/\sqrt{27}$, because all $k$ missing powers divide the predecessor and are assigned to $C'$. Together with $U\le K H^{nk}$ this requires $P<K^{3/(2k)}H^{3n/2}$, contrary to the first threshold.

At the product, all three parts vanish at all $m$ assigned nodes. Using the largest assigned node as origin, each part has $m-1$ positive roots, and the total also has those roots. There must be at least $3\lceil(m-1)/2\rceil$ distinct negative positions and at least $1+\lfloor(m-1)/2\rfloor$ positive positions of the total, all in degrees 1 through $n$. The preceding digit comparison makes these two sets disjoint. This gives $m\le b$ and $2n-3b>0$.

Each full power in $S$ divides at least two of the $L_i(P)$, since it divides two actual parts and is coprime to $\Phi(P)$. Therefore $S^2\le K^3 H^{3nk}$. Also $U\le d(2P)^m$ and $N-1>dP^n$. Applying $27T^2\le U^3$ gives
$$27P^{2n-3m}< 2^{3m}d K^3 H^{3n(n-m)}.$$
The exponent $3n(n-m)/(2n-3m)$ increases with $m$ on this range, and $C\ge 2^{3m}d K^3$. This proves the second threshold.

For every fixed $n,d$ both ranges contain infinitely many choices of distinct primes. Partition $[Y,2Y]$ into intervals of length a sufficiently small constant times $Y^{1/\gamma}$ for the respective exponent $\gamma$. The prime number theorem and convexity give $\gg_{n,d}Y^{1+(n-1)/\gamma}/(\log Y)^n$ choices of $n$ primes in the same interval: bins with fewer than half the mean occupancy contain at most half the total prime count, and one can apply the power-mean inequality to the others. The primes in each choice are recovered from $F(P)$ by unique factorization.
