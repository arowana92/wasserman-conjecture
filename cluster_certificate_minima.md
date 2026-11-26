# How often both weights need all four primes

There is an absolute constant $c>0$ such that, for large $X$ and every integer $1\le D\le c\log X$, at least a constant times
$$\frac{X^{2/3}}{D(\log X)^3}$$
distinct products $M=PQR\le X$ have the following property. For **every** $1\le d\le D$, one list of four primes certifies all choices of three positive trinomial rows at both $dM$ and $dM+1$. Four is the minimum size at each of the two weights. The primes $P,Q,R$ belong to every list for that product. The fourth prime can depend on $d$.

There are two versions for longer ranges. Fix $\epsilon>0$. In each version the same lower bound on the number of products holds, with an implied constant allowed to depend on $\epsilon$.

* For $D\le X^{31/1269-\epsilon}$ there is a shared four-prime list for **every** $d\le D$. For each product both minimum sizes equal four at all but $O_\epsilon(D/\log X)$ multipliers.
* For $D\le X^{31/669-\epsilon}$ there is a shared four-prime list of minimum size at both weights for all but $O_\epsilon(D/\log X)$ multipliers of each product. The [three-prime theorem](clustered_predecessors.md) still proves the conjecture at both weights for every multiplier of these products.

In either version we can simultaneously require the shared list to be optimal at every $d\le\min(D,c_\epsilon\log X)$ for a suitable $c_\epsilon>0$.

## Prime forms that could shorten the list

Use $Y=X^{1/3}/2$, $H=\lfloor\sqrt{Y/D}/20\rfloor$, and consider all prime triples $Y\le P<Q<R\le2Y$ with $R-P\le H$ and $R-Q>Q-P$. Their number is $\gg B$, where $B=YH^2/(\log Y)^3$, by the [short-interval count](clustered_predecessors.md#counting-whole-intervals-of-multipliers). In the two power ranges above $\log H\asymp\log Y$ uniformly. Constants for the logarithmic range can be chosen absolute.

At a weight having a four-prime list, its minimum size is one or four. In the absence of a universal prime, for each prime in a list of at most three one can select a row avoiding it. A prime $b$ is universal exactly when the base-$b$ digit sum of the weight is at most two. If $b>\sqrt N$ at weight $N$, the only possibilities for large $N$ are $N=b,2b,b+1$. The number of weights $N\le T$ with any universal prime $b\le\sqrt N$ is $O(\sqrt T\log^2 T)$ by counting bases and the positions of the one or two digit units.

Write $n=dPQR$. The three large prime factors rule out $n$ being prime or twice a prime. The only remaining large-base possibilities at the two weights are primality of one of
$$n-1,\qquad n+1,\qquad (n+1)/2,$$
where the last is considered only if integral. For each fixed $d$, a weight specifies at most two products when the choice between $n$ and $n+1$ is unspecified. Thus the number of incidences $(M,d)$, $d\le D$, involving a small universal base is $O(Y^{3/2}D^{3/2}\log^2(DY))$. Throughout both power ranges this is $o(B/\log Y)$: it suffices that $D^5$ be smaller than $Y$ by a fixed power, and $3\cdot31/669<1/5$.

## A uniform upper sieve for the large bases

Fix $Q,R$ and an integer $d$. In each of the three cases above, no odd prime below $Y$ divides $P$ or the relevant form $dPQR\pm1$, since the prime being tested exceeds $Y$. Sieve the possible $P$ in an interval of $H$ consecutive integers using odd primes $\ell\le z$ with $z<Y$. Remove the classes
$$P\equiv0\pmod\ell,\qquad dQRP\equiv\mp1\pmod\ell.$$
There are $\nu(\ell)=2$ forbidden classes if $\ell\nmid d$ and one otherwise. Omitting 2 accommodates the form that is twice a prime. The upper bound for either sign is
$$\ll \frac{H}{(\log H)^2}\, f(d),\qquad
 f(d)=\prod_{\substack{\ell\mid d\\\ell>2\ {\rm prime}}}\frac{\ell-1}{\ell-2}. \tag{1}$$

Here are details giving uniformity in $d,Q,R$. Use the Selberg upper bound with $A=H^{1/4}$, $z=A^{1/16}$ and squarefree moduli supported on the odd primes through $z$. Write $g(\ell)=\nu(\ell)/\ell$ and $h(\ell)=g(\ell)/(1-g(\ell))$, and extend multiplicatively to those moduli. The bound is $H/J$ plus an error $O(H^{1/2+o(1)})$, where $J=\sum_{m\le A}h(m)$ over these moduli. Indeed each squarefree $m$ has at most $2^{\omega(m)}$ forbidden CRT classes in the intersection, so the interval remainder has that same bound. The error sum has $m\le A^2$, with coefficient bounded by the ternary divisor function. [Tao, *Some sieve theory*, Theorem 30](https://terrytao.wordpress.com/2015/01/21/254a-notes-4-some-sieve-theory/)

Normalize the weights $h(m)$ on all such squarefree moduli to probability one. Each prime occurs independently with probability $g(\ell)$. Hence the expected $\log m$ is at most $2\sum_{\ell\le z}(\log\ell)/\ell\le4\log z$: use $\lfloor z/\ell\rfloor\ge z/(2\ell)$ and $\sum_{\ell\le z}\lfloor z/\ell\rfloor\log\ell\le\log(\lfloor z\rfloor!)$. Markov's inequality gives $J\ge(3/4)\prod_{\ell}(1-g(\ell))^{-1}$. For $\ell>2$, $1-2/\ell\le(1-1/\ell)^2$. Replacing two forbidden classes by one costs $(\ell-1)/(\ell-2)$. Finally $\prod_{\ell\le z}(1-1/\ell)^{-1}\ge\sum_{j\le z}1/j$. This proves (1).

There are $O(YH/(\log Y)^2)$ choices of $Q,R$ by the prime number theorem and the short-interval Brun–Titchmarsh bound. For each such pair $P$ lies in an interval of at most $H$ integers. Consequently the large-base incidences at each $d$ number $O(B f(d)/\log Y)$. Expanding the product in (1) gives
$$\sum_{d\le E} f(d)\le E\prod_{\ell>2}\left(1+\frac1{\ell(\ell-2)}\right)\ll E.$$
The product converges by comparison with $\sum 1/j^2$. Thus the incidence count through any $E\le D$ is $O(BE/\log Y)$, apart from the small-base incidences already bounded.

## Combining the sieve with the fourth-prime count

For a given $d$, call a product a failure here when the interval $[dM-Y/2,dM]$ contains no prime. The [close-product and prime-gap estimates](shared_cluster_certificates.md#finding-the-fourth-prime-at-all-multipliers) bound failures by
$$\ll_\eta (Y^{1+\eta}H^2)^{1/2}\,Y^{-1}(dY^3)^{123/200+\eta}
\ll_\eta Y^{369/200+O(\eta)}D^{-1/2}d^{123/200+\eta}.$$
This uses the power $123/100$ in the mean square gap estimate and Cauchy–Schwarz on bins. At every other product a prime in that interval completes the shared list.

Summing through $D$ gives $O_\eta(Y^{369/200+O(\eta)}D^{223/200+O(\eta)})$ incidences. In the smaller power range this is $o(B)$ because $D^{423}<Y^{31}$ with a power saving. In the larger range it is $o(BD/\log Y)$ because $D^{223}<Y^{31}$ with a power saving. Take $\eta$ sufficiently small in each case. Summing only through $E\le c_\epsilon\log X$ gives $o(B)$ even in the larger range: the ratio to $B$ is at most $Y^{-31/200+O(\eta)}\sqrt D$ times a power of $\log Y$.

Start with at least $c_0B$ products. Discard products with any failure or universal prime through $\min(D,c_\epsilon\log X)$. The early prime-gap and small-base costs are $o(B)$, and taking $c_\epsilon$ sufficiently small makes the early large-base cost less than $c_0B/4$. In the smaller power range we can discard all products with a prime-gap failure anywhere, costing another $o(B)$. The remaining incidences of exceptional multipliers total $O_\epsilon(BD/\log X)$ in either version. Discarding products with more than a sufficiently large constant times $D/\log X$ such multipliers leaves $\gg_\epsilon B$ products. This proves both power-range statements and, on taking $D\le c\log X$ with one fixed $\epsilon$ in the smaller range, the absolute-constant theorem.

## An explicit initial interval

Take $P=10^{46}+121$, $Q=10^{46}+537$ and $R=10^{46}+117073$. Write $M=PQR$. At **every one of the 1,123 multipliers $d=1,\dots,1123$** the two weights $dM,dM+1$ have a shared list of four primes and both minimum sizes are four. The largest distance from $dM$ to the chosen fourth prime is 2,876, at $d=213$. At $d=1124$ the run ends because $dM-1$ is prime and alone certifies $dM$.

Run `python3 check_cluster_minima.py` with python-flint installed. It verifies the core primes, finds and proves a fourth prime below $dM$ for each multiplier, and checks minimum size in two ways. The independent check allows composite bases in its digit search, uses elementary integer bisection for every possible leading exponent, and exhibits a nontrivial gcd or a Fermat compositeness witness for each candidate base. Small weights check this enumeration by direct digit division.
