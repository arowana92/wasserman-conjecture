# Infinite families with bounded prime certificates

**Theorem.** For all sufficiently large $X$, at least $cX/(\log X)^2$ distinct odd composite integers $N$ in $[X,2X]$ have the following two properties, for an absolute $c>0$.

1. Some set of at most seven primes certifies Wasserman’s assertion for all three-row selections at $N$.
2. No prime divides every positive trinomial coefficient at $N$.

The primes in the certificate depend on $N$. The proof uses the Bombieri–Vinogradov theorem, the lower-bound linear sieve and Guth–Maynard’s almost-all short-interval theorem. Only the frequency estimate uses analytic number theory.

## A sufficient condition that can be checked for one weight

Let $N=2q+1$ with $q$ an odd prime. For each prime $p$ dividing $N$, let $P_p=p^{v_p(N)}$. Suppose that $r$ is a prime with $q<r\le N$ and
$$N-r<2\min_{p\mid N}P_p.$$
Then the distinct primes in $\{q,r\}\cup\{p:p\mid N\}$ certify $N$.

Indeed, suppose three rows avoid every candidate prime between them. A row avoiding $q$ must assign the three base-$q$ digit units of $N$ to three nonempty summands. It is $(1,q,q)$ up to permutation. A row avoiding $r$ must have a part at least $r$, since $r\le N<2r$. It is therefore a different row. Its two smaller positive parts sum to at most $N-r$. For any $p\mid N$, a row avoiding $p$ has all parts divisible by $P_p$, because $N$ ends in $v_p(N)$ base-$p$ zeroes. Neither of the first two rows can do so, one because it contains 1 and the other because its two smaller parts would sum to at least $2P_p$. Thus the third row must avoid every $p\mid N$. Each of its parts would be a positive multiple of $N$, a contradiction. The argument also covers repeated rows.

For example, $N=119=7\cdot17=2\cdot59+1$ and $r=113$ give the four-prime certificate $\{7,17,59,113\}$.

## Constructing many eligible composite weights

Throughout, the decimal exponents denote exact rational numbers. Put
$$a=0.141,\quad b=0.145,\quad \gamma=0.176,\quad \delta=0.354,\quad h=0.137.$$
Take $X$ large, $z=X^\gamma$, $D=X^\delta$. Let $\mathcal M$ consist of primes $m$ in $[X^a,X^b]$. Form a multiset of integers $u=(2q+1)/m$ over all pairs of primes $(q,m)$ with $X<2q+1\le2X$, $m\in\mathcal M$, and $m\mid 2q+1$. We claim that $\gg X/(\log X)^2$ members, counted with multiplicity, have no prime factor below $z$.

**Distribution estimate.** Let $T$ be the logarithmic-integral mass of the interval $((X-1)/2,(2X-1)/2]$ and put $Y=T\sum_{m\in\mathcal M}1/(m-1)$. By the prime number theorem and partial summation the reciprocal sum tends to $\log(b/a)>0$, so $Y\asymp X/\log X$. For squarefree $d$ let $A_d$ count members of the multiset divisible by $d$, and let $\rho(d)=1/\varphi(d)$ for odd $d$ and zero for even $d$. We have, for every fixed $B>0$,
$$\sum_{d\le D}\mu^2(d)|A_d-Y\rho(d)|\ll_B X/(\log X)^B.$$

For odd $d$, the condition for a given $m$ is $2q\equiv-1\pmod{md}$, one reduced residue class even when $m\mid d$. The modulus is at most $X^{b+\delta}=X^{0.499}$. In the double sum over $m,d$ each modulus occurs a bounded number of times, since each possible distinguished factor $m$ is prime and at least $X^a$. Applying Bombieri–Vinogradov in its prime-counting form at the two endpoints replaces each progression count by $T/\varphi(md)$ with the stated total absolute error. For $m\nmid d$ the main term is $T/((m-1)\varphi(d))$. If $m\mid d$ the absolute difference is $T/(m(m-1)\varphi(d))$. Summing this correction even without the restriction $m\mid d$ costs
$$\ll T\sum_{m\ge X^a}\frac1{m(m-1)}\sum_{d\le D}\frac1{\varphi(d)}
\ll TX^{-a}\log(2D).$$
Here the last sum is $O(\log(2D))$, by summing the identity $d/\varphi(d)=\sum_{v\mid d}\mu^2(v)/\varphi(v)$ and the convergent series $\sum_v\mu^2(v)/(v\varphi(v))$. For even $d$ both counts are zero.

**Lower-bound sieve.** The local densities are $\rho(2)=0$ and $\rho(p)=1/(p-1)$ at odd primes. They have sieve dimension one by Mertens’ product formula, since the ratio of $1-1/(p-1)$ to $1-1/p$ equals $1-1/(p-1)^2$. The sieve product up to $z$ is consequently $\asymp 1/\log z$. Our nonnegative multiset is supported on $u\le X$ for large $X$ and has the distribution estimate above. The lower linear sieve applies with $s=\log D/\log z=0.354/0.176>2$, so its lower function is strictly positive. It gives the claimed $\gg Y/\log z$ sifted members. A statement with precisely these distribution hypotheses is Definition 2.4 and Lemma 2.5 of *Weighted sieves with switching*. [Linear sieve and its hypotheses](https://www.cambridge.org/core/journals/mathematical-proceedings-of-the-cambridge-philosophical-society/article/weighted-sieves-with-switching/986429394BA969687D48224E177E2F1C)

For every surviving pair, $u>1$ and $m<z$, hence $m\nmid u$. Also $u$ has at most four prime factors counted with multiplicity, since $m z^5\ge X^{a+5\gamma}=X^{1.021}>2X$ for large $X$. Thus $N=mu$ has at most five prime factors and has the distinguished factor $m$ exactly once. Every prime factor of $N$ is at least $X^a$. The representation by a surviving pair is unique, because $m$ is the only prime factor below $z$ and $q=(N-1)/2$. This gives $\gg X/(\log X)^2$ distinct composite weights.

**A prime near each retained weight.** Put $H=\lfloor X^h\rfloor$. Guth–Maynard’s Corollary 1.4 gives a prime in $(x,x+H]$ for all but $O(X\exp(-c_1(\log X)^{1/4}))$ integer starts $x\in[X,2X]$, since $h>2/15$. Take $x=N-H$. The lower boundary costs at most $H+1$ weights, and the map from $N$ to $x$ is injective. Both losses are $o(X/(\log X)^2)$. For every retained weight the prime $r$ satisfies $r>q$ and $N-r<H<X^a$ for large $X$. The elementary criterion now gives a certificate of size at most $2+5=7$. [Guth–Maynard, Corollary 1.4](https://ora.ox.ac.uk/objects/uuid%3Aad11b8bf-ad2b-4ebf-a627-647f023c378f/files/s73666708v)

Finally, a universal prime would have digit sum at most two for $N$ in its base. An odd $N$ of digit sum two in an odd base cannot occur, since it would be a sum of two odd prime powers, allowing exponent zero. A digit sum of one would make $N$ a prime power, which is excluded by $m\parallel N$ and $N>m$. In base two, odd $N$ of digit sum at most two and $N>1$ has $N-1$ a power of two. Here $N-1=2q$ with $q$ an odd prime. This completes the proof.
