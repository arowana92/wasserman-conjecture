# A single prime list across weights needs at least square-root order up to logarithms

Let $\mathcal P$ be any list certifying every weight $N=2p\le X$, $p$ prime. Each of these weights individually has a one-prime certificate. Nevertheless, for all sufficiently large $X$,
$$
 \#\{r\in\mathcal P:r>(X/2)^{3/11}\}
                            \gg \frac{\sqrt X}{\log X}.                 \tag{1}
$$
There are even $\gg\sqrt X/\log X$ different weights with explicit selections whose common gcds are pairwise coprime and supported on primes above $(X/2)^{3/11}$. At a higher prime cutoff $c\sqrt X/\log X$, there are still $\gg\sqrt X/(\log X)^2$ such weights for an absolute $c>0$.

A finite construction gives **45,287** different weights between five and ten million with pairwise disjoint sets of possible common primes for the displayed selections. All their common primes exceed 2,236. Thus a single list certifying these weights must contain at least 45,287 primes above 2,236.

## The two possible common primes

At $N=2p$ every positive trinomial coefficient contains $p$: there are only two digit units in base $p$ for odd $p$, so a positive three-way split must carry. The case $p=2$ is immediate as well. Suppose $2p-1=M$ has at most two prime factors counted with multiplicity and $s_2(N)\ge3$. Take the endpoint and an odd coefficient, the latter obtained by allocating the binary digit units of $N$ among three nonempty parts. Any common prime is now either $p$ or a prime of $M$.

If $M=rq$ with $r\ne q$ primes, choose the third row to avoid $r$. This is possible because $s_r(N)=1+s_r(q)\ge3$: the distinct prime $q$ is not a power of $r$. The common primes of the three rows then belong to $\{p,q\}$. Exchanging $q$ and $r$ gives another triple with common primes in $\{p,r\}$. Each triple has $p$ as an actual common divisor. If $M=q$ or $q^2$, the endpoint and the odd row already leave only $\{p,q\}$, and the odd row may be repeated.

Represent each possible two-prime set by an edge between its primes. Every certifying list must hit every edge. When $M=rq$ is squarefree with two factors, omitting $p$ requires including **both** $r$ and $q$.

## Counting the weights

Use the Chen–Friedlander–Iwaniec theorem for linear forms $ap+b=cs$ with $(a,b,c)=(2,-1,1)$. These coefficients are pairwise coprime, $a,c\ge1$, $b\ne0$, and $2\mid abc$. For $t=X/2$ it gives $\gg t/(\log t)^2$ primes $p\le t$ such that $2p-1$ has at most two prime factors, each exceeding $t^{3/11}$. The theorem permits negative $b$. This is the unconditional result stated as Theorem 6 in the cited paper. [Garcia, Luca, Shi and Udell, Theorem 6](https://arxiv.org/pdf/1906.05927)

Discard the $O((\log t)^2)$ weights with at most two binary digit units and the finitely many trivial inputs. Also discard the $O(t^{3/11})$ inputs with $p\le t^{3/11}$. Every vertex of the remaining graph now exceeds $t^{3/11}$.

If the edges are hit by $k$ vertices, at most $k$ of the eligible weights have their own $p$ among the chosen vertices. Among the others with $M=q$ or $M=q^2$, the vertex $q$ must be chosen and determines at most two weights. For $M=rq$ with distinct factors the unordered pair $\{r,q\}$ must be chosen and determines $p$ uniquely. The total count is therefore at most $3k+\binom{k}{2}$. Comparison proves (1). The endpoints of a maximal matching hit every edge, so that matching has at least half the lower bound. Its edges come from different weights because the edges supplied by one weight all contain its $p$. The corresponding actual coefficient gcds are pairwise coprime.

For the higher cutoff discard also $p<t/(\log t)^2$ at cost $O(t/(\log t)^3)$ by the prime number theorem. Keep only the edge $\{p,q\}$ where $q$ is the largest prime factor of $M$. Then $q\ge\sqrt M\gg\sqrt t/\log t$. Any vertex accounts for at most one weight as the weight prime $p$. As the largest factor $q$, it accounts for at most one with $M=q$, and the others have $M=rq$ for a prime $r\le2t/q\ll\sqrt t\log t$. There are $O(\sqrt t)$ such primes, again by the prime number theorem. Removing all edges incident to the endpoints of one edge costs $O(\sqrt t)$ edges. A greedy matching therefore has the claimed second size. Repeated unordered edges, if any, can be counted with their weights in this argument. [Selberg, prime number theorem](https://www.math.lsu.edu/~mahlburg/teaching/handouts/2014-7230/Selberg-ElemPNT1949.pdf)

## Checking finite witnesses

The [checker](check_double_prime_pool.py) enumerates primes with $X/2<2p\le X$ by a sieve. It checks the factorization of $2p-1$, constructs an odd row and, when necessary, an $r$-avoiding row, and evaluates their full gcd with the endpoint by Legendre's factorial floors. In decreasing order of $p$ it keeps the edge to the largest prime factor of $2p-1$ when both vertices are unused. This matching is a rigorous lower bound, without a claim of optimality.

| $X$ | Eligible weights in the upper half | Pairwise disjoint edges |
|---:|---:|---:|
| 4,000 | 68 | 63 |
| 10,000,000 | 54,166 | 45,287 |

The small run also multiplies the exact binomial coefficients for every retained edge and compares the integer gcd to the factorial-floor result.

~~~sh
python3 check_double_prime_pool.py
~~~
