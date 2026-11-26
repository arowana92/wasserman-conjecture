# A prime divisor certified by the scale of a row

Let $N=VU$, with $V,U$ positive integers. Consider every positive trinomial row of sum $N$ whose parts are multiples of $V$. Let $w=\ell^e>1$ divide $N-1$, with $\ell$ prime. For any $m=\ell^j$, $j\ge1$, write $[x]_m$ for the least nonnegative residue. A necessary condition for such a row to avoid $\ell$ is
$$[-V]_m\left\lfloor\frac{U-1}{w}\right\rfloor
\ \ge\ 2m-\left[\left\lfloor\frac Nw\right\rfloor\right]_m. \tag{1}$$
Failure at any one block certifies $\ell$ as a divisor of **every** coefficient in this class. In particular if $V\equiv-1\pmod w$ and $U\le w^2$, every coefficient is divisible by $\ell$.

Indeed a carry-free row has residues $0,0,1$ modulo $w$. Since $V$ is a unit, the two zero parts are $Vw x,Vw y$ for positive integers $x,y$ with $x+y\le\lfloor(U-1)/w\rfloor$. Their next $j$ base-$\ell$ digits have values $[Vx]_m,[Vy]_m$. Their sum is at most the corresponding block of $N$, by carry freedom. Put $J=[-V]_m>0$. Then $[Vx]_m\ge m-Jx$ for every $Jx>0$, including $Jx\ge m$. Adding the corresponding inequality for $y$ proves (1). For $V\equiv-1\pmod w$ take $m=w$. The right side is at least $w+1$.

This applies even before individual rows or CRT assignments are enumerated. For example, let $P,Q$ be full powers of distinct primes and let $w$ be a prime power dividing $\gcd(P+1,Q+1)$. If $\max(P,Q)\le w^2$, the full conjecture holds at $N=PQ$. By the [endpoint theorem](all_weights_endpoint_reduction.md) it suffices to treat selections containing $(1,1,N-2)$. In a possible counterexample the two other rows would have to avoid the base primes of $P,Q$ separately, since a row avoiding both would have all parts divisible by $N$. Each of those rows has parts divisible by its assigned prime power, so (1) forces the base of $w$ in both coefficients. It divides the endpoint coefficient as well.

## Infinitely many certified products of two primes

For every $\epsilon>0$ and all sufficiently large $Y$, at least $\gg_\epsilon Y^{3/4-\epsilon}$ distinct weights $N\le Y$ are products of two odd primes $P,Q$ with a prime $w\mid\gcd(P+1,Q+1)$ and $\max(P,Q)<w^2$. For each of them the three primes $P,Q,w$ certify every pair of rows together with the endpoint. The endpoint theorem then proves the full conjecture at all these weights.

Fix $1/2<\theta<0.51$. Use the lower bound in the average Brun–Titchmarsh theorem of Baker and Harman at the fixed residue $-1$. For all but $O_\theta(x^\theta/\log^2 x)$ moduli $w$ in $(x^\theta,2x^\theta]$, restricted here to primes, there are $\gg_\theta x^{1-\theta}/\log x$ primes $P\le2x$ in this residue class. One can read this directly from the minorant formulation. Its total mass is $\gg_\theta x/\log x$, and restricting to integers coprime to a prime $w$ removes at most one in net mass because the minorant is bounded above by the prime indicator. Its sum of absolute progression errors is $O_{A,\theta}(x/\log^A x)$. Comparing each error with a small constant times $x^{1-\theta}/\log x$ gives the claimed exceptional set. This use of the published distribution estimate concerns moduli beyond $\sqrt x$. [Baker and Harman, *The Brun–Titchmarsh theorem on average*](https://doi.org/10.1007/978-1-4612-4086-0_4) and [minorant formulation, introduction and Theorem 2.1](https://arxiv.org/html/2602.20917v1#S2).

The prime number theorem supplies $\gg_\theta x^\theta/\log x$ good prime moduli. At each, discard the at most three primes $P\le3w$ in the progression and take all unordered pairs of remaining primes. Both normalized sums now pass the elementary size test $2w+1\le\min(P,Q)$. For large $x$ the product of any two moduli in the range exceeds $2x+1$, so no prime $P$ can occur for two different moduli. Unique factorization therefore gives $\gg_\theta x^{2-\theta}/\log^3 x$ distinct weights at most $4x^2$. Take $x=\sqrt Y/2$ and choose $\theta$ sufficiently close to $1/2$.

Three is the minimum size for the endpoint-pair certificate at these weights. Each prime $r$ is individually avoidable by splitting at least three base-$r$ digit units. For $r=P$ the digit sum of $N$ is that of the distinct odd prime $Q$ in the odd base $P$. It is odd and cannot be one. The argument at $Q$ is the same. At a prime not dividing $N$, a digit sum of at most two would require $N-1$ to be a prime power, whereas both $2$ and the odd prime $w$ divide it. Thus two rows suffice to avoid any prescribed pair of primes between them.

## Four primes for an 81-digit weight

Put $w=2^{127}-1$, $P=80w-1$, $Q=282w-1$, and $N=PQ$. The numbers $P,Q,w,R=N-94$ are prime. These four primes form a minimum certificate for the **original three-row problem** at this weight. This certificate also gives a direct proof at $N$, independent of the general endpoint reduction.

No positive row can avoid two of the four primes:

- Avoiding $P$ and $Q$ together would require three positive multiples of $N$.
- Avoiding $P$ or $Q$ forces the corresponding scale. Both are $-1\bmod w$ and the complementary sums are below $w^2$, so (1) forces $w$.
- A row avoiding $R>N/2$ has one part at least $R$ and its two remaining positive parts sum to at most 94. They cannot both be multiples of $P$ or $Q$. If the row avoided $w$, at least one of them would have to be a multiple of $w>94$.

Any three rows can therefore avoid at most three primes from the set. For minimality, every individual prime $r$ is avoidable at $N$. It suffices that the sum of the base-$r$ digits of $N$ is at least three, by distributing digit units to three nonempty parts. At $r=P,Q$ this holds since $3\le\lfloor Q/P\rfloor<P$, $Q<P^2$, and $P<Q$. At $r\nmid N$ a digit sum at most two for $N>2$ would require $N=1+r^i$, since its units digit is nonzero. But $N-1$ is divisible by both 2 and $w$. Thus from any set of at most three primes one can choose a row avoiding each of its primes.

## Two square-root restrictions in the close-gap families

Use the notation of the [six residual families](close_power_gaps.md). In a putative counterexample put $s=p^b$, $t=p^v$, $G=p^{h-b}+1$, $L=(k+1)s$, $D=LG+t$. If $t\mid s$, **each full prime power $w$ in $D$ satisfies $w^2<D/t$**. Indeed the companion scale is $V=N/G$ and $tV\equiv-L\pmod w$, so $[-V]_w\le L/t$. With normalized sum $G$, (1) gives $w+1\le(L/t)\lfloor(G-1)/w\rfloor$, which implies the assertion.

If $t\mid ks$, **each full prime power $w$ in $G$ satisfies $w^2\le LG$**. Suppose $w^2>LG$, and write $g=G/w$. If the base prime of $w$ is assigned to the first row, its two unshifted parts divided by $w$ have residues $kz,z\bmod w$, where $tz\equiv-sg\bmod w$ and $0\le z<w$. Since $w>L g$, one has $tz\ge w-sg$ and $[kz]_w=w-ksg/t$. The residues sum to more than $w$, forcing a carry. Thus that prime is absent from the first row's assigned factor $U\mid G$, whence $U\mid g$. But $g<\sqrt{G/L}$ gives $U^3<27D^2$, contradicting the required size bound. In particular, in each of the three residual families with $v=0$, every full prime power of **both** $G$ and $D$ is at most $\sqrt D$.

With python-flint installed, run `python3 check_scaled_row_blocks.py`. It checks (1) against exact binomial products in 289,365 comparisons through weight 150, including 110 obstructions that pass the test $U\ge2w+1$. It also verifies all four large primes by recursive modular order certificates. For each tested integer $n>2$, the code proves the factorization of $n-1$ recursively. For each prime $q\mid n-1$, its witness has order divisible by the entire $q$-part of $n-1$ modulo every prime divisor of $n$. Thus every such divisor is $1\bmod(n-1)$, proving primality.
