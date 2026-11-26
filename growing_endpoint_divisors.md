# A growing common divisor for every pair of trinomials

Outside $O(X^{0.9993})$ weights $N\le X$, **every pair** of positive trinomial rows $A,B$ satisfies at least one of
$$
 \gcd(N(N-1),K_A,K_B)>N^{1/8},\qquad
 2^{\min(v_2(K_A),v_2(K_B))}>N^{1/8}, \tag{1}
$$
where $K_A,K_B$ are their coefficients. In particular the part of their gcd supported on primes of $N(N-1)$ exceeds $N^{1/8}$, allowing the exponents of these primes in the gcd. If one coefficient is odd, the **actual gcd with the endpoint coefficient** exceeds $N^{1/8}$. There are no restrictions on the sizes or repetitions of the parts.

There is also a statement about individual common primes. For $X\le N\le 2X$ and $2\le y\le X$, outside
$$ O\!\left(X^{0.9993}+X\exp\!\left[-\frac{\log X}{32\log(2y)}\right]\right) \tag{2} $$
weights, **every** pair with $2^{\min(v_2(K_A),v_2(K_B))}\le N^{1/8}$ shares a prime $p>y$ dividing $N(N-1)$. Thus for any $y=X^{o(1)}$ growing to infinity, this holds on a set of relative density tending to one, including all pairs with one binary carry-free row. For triples containing $(1,1,N-2)$, that same prime divides all three coefficients. The general three-row conjecture remains open.

## Counting after removing common prime powers

Write $M=N-1$, $D=\gcd(N,K_A,K_B)$ and $J=\gcd(M,K_A,K_B)$. Use the [partial-power assignment](endpoint_pair_bound.md#retaining-partial-prime-powers), splitting $N=DUV$, $M=JST$. Each retained prime power goes to a row of smaller valuation at its prime base. Name a row of smaller binary valuation $A$; choose it for 2 in a tie. Its parts are multiples of $U$, and modulo each retained power in $S$ its three residues are $0,0,1$. The corresponding factors of $B$ are $V,T$. The removed factors may share primes with the factors retained on the same side of $N,M$.

On $Q/2\le N<Q=2^L$, for $1\le Y,Z\le Q$ the number of weights permitting an assignment with $US\le Y$, $DJ\le Z$ is
$$ \ll Q^{-2/5+o(1)}Y^{8/5}Z^{13/5}. \tag{3} $$
Here $Q^{o(1)}$ is uniform. Adapt the [square-modulus count](common_factor_count.md#counting-small-assigned-moduli-using-one-predecessor-factor), fixing $D,J$. For the smallest of the three bins in $T=dt$ the normalized parts of $B$ at the designated position and the other two positions combined are $zt,dy$, respectively. Thus
$$ DU=zt+dy,\quad V(zt+dy)=JSdt+1,\quad Vzt\equiv1\pmod d,\quad d\le T^{1/3}. $$
In a dyadic box write $U\asymp U_0$, $S\asymp S_0$, $d\asymp d_0$, $K_0=U_0 S_0$, $V_0=Q/(DU_0)$ and $T_0=Q/(JS_0)$. Given $S,V,d,z$, the weight lies in one class modulo $VJSd^2$ since $\gcd(V,JSd)=1$. The choice of $z$ costs $O(DU_0d_0/T_0)$. This gives the box count
$$ \ll Q^{2/3}J^{1/3}S_0^{4/3}+DK_0. $$
Instead fix $S,d,z,y$, with $y\ll DU_0/d_0$. The positive integer $JSd^2y-z=(JSd-Vz)DU$ determines $DU$ up to a divisor choice and then determines $t,N$. Positivity follows from $(JSd-Vz)t=Vdy-1>0$: $Vdy$ is the sum of two positive actual parts, including when $V=1$. The second count is
$$ \ll Q^{-2/3+o(1)}D^2 J^{2/3}K_0^2 S_0^{-1/3}. $$
Taking powers $1/5,4/5$ as in the square-modulus count gives $Q^{-2/5+o(1)}D^{8/5}J^{3/5}K_0^{8/5}$, plus $DK_0$. The latter is absorbed since $27T^2\le(DU)^3$ implies $DJK_0\gg Q^{2/3}$. Sum over the boxes and at most $O(Z(1+\log Z))$ pairs $D,J$ with $D^{8/5}J^{3/5}\le Z^{8/5}$. This proves (3).

## A binary digit count that allows carries

To prove (1), count weights with $DJ\le Q^{1/8}$ and $v_2(K_A)\le L/8$. Write $s_2$ for binary digit sum. Legendre's formula gives the exact identity
$$ \sum_i s_2(A_i)=s_2(N)+v_2(K_A). \tag{4} $$
This controls the entropy cost of the three parts even when their digit supports overlap. Use $r=s_2(N)/L$, fix the positions of parts equal to one, and choose a canonical row pair and assignment. For $r>517/1000$ the binomial upper tail gives at most $2^{L H(517/1000)}$ weights, with $H$ the binary entropy.

Suppose first all three parts of $A$ exceed one. For $US\le Q^{33/50}$ use (3). Otherwise, on the uniform distribution on the chosen weights, the label $P=(U,S)$ leaves at most $O(Q^{1-33/50})$ choices for $N$ by the ordinary CRT. Each part $x\ge2$ of $A$ specifies $P$ up to $Q^{o(1)}$ possibilities, from $U\mid x$ and $S\mid x(x-1)$. Apply the [shared-label entropy inequality](common_factor_count.md#information-that-every-individual-part-supplies) with three parts. By (4), subadditivity over digits, concavity of $H$, and its increase below $1/2$, the mean of the three part entropies is at most $L H((517/1000+1/8)/3)$. The exponent for the count is at most
$$ H((517/1000+1/8)/3)+(2/3)(1-33/50)+o(1). $$

If a part of $A$ equals one, $U=1$ and every retained power in $S$ divides both other parts because the exceptional residue is already occupied. For the other row $27T^2\le D^3$, so
$$ S=M/(JT)\ge 3^{3/2}(Q/2-1)/(J D^{3/2})\ge Q^{1-3/16}/2 $$
for all sufficiently large $Q$ (also directly using $T\le(D/3)^{3/2}$). In particular a second part equal to one is impossible. Apply the same label inequality to the two other parts, whose sum determines $N$ by adding one. Each identifies $S$ up to divisor choices. The conditional entropy of $N$ given $S$ is at most $(3/16)L+O(1)$, whereas the mean of their entropies is at most $L H((517/1000+1/8)/2)$ by (4). The exponent in this branch is at most
$$ H((517/1000+1/8)/2)+3/32+o(1). $$
Both displayed exponents, the binomial tail exponent, and the exponent $-2/5+(8/5)(33/50)+(13/5)/8$ in (3) are strictly below $0.9993$. This proves (1) on summing dyadic blocks.

## From a common divisor to a common prime

Let $m_y(n)$ be the product of full powers in $n$ of primes $\le y$. Uniformly in $y\ge2$, with $\sigma=1/(2\log(2y))$, we have
$$ \sum_{n\le x} m_y(n)^\sigma\ll x. $$
Indeed expand $m_y(n)^\sigma$ as a divisor sum with nonnegative multiplicative coefficients supported on $y$-smooth integers. At $p^j$, $p\le y$, the coefficient is $p^{j\sigma}-p^{(j-1)\sigma}$. Counting multiples of each divisor and extending the sum gives the upper mean
$$ \prod_{p\le y}\left(1+\frac{p^\sigma-1}{p-p^\sigma}\right). $$
Its logarithm is $O(\sigma\sum_{p\le y}(\log p)/p)=O(1)$. Indeed $p^\sigma\le\sqrt e<2$ and the numerator is $O(\sigma\log p)$. For $m=\lceil y\rceil$, the elementary bound $\sum_{p\le m}(\log p)/p\le2\log m$ follows by $\lfloor m/p\rfloor\ge m/(2p)$ and $\sum_{p\le m}\lfloor m/p\rfloor\log p\le\log(m!)$.

If the first gcd in (1) exceeds $N^{1/8}$ and all its prime factors are at most $y$, at least one of $m_y(N),m_y(N-1)$ exceeds $X^{1/16}$ in $[X,2X]$. The preceding moment bounds the number of such weights by $O(X\exp[-\log X/(32\log(2y))])$. This proves (2).

Run the standard-library verifier to check the strict margins, the partial-power congruence on exact multinomial pairs and the divisor-sum identity.

~~~sh
python3 check_growing_endpoint_divisors.py
~~~
