# A predecessor prime dividing one trinomial coefficient

**Finite theorem.** For every $3\le N\le2{,}000{,}000{,}000$, every odd positive trinomial coefficient $K=N!/(a!b!c!)$ has $\gcd(K,N-1)>1$. Equivalently, every positive trinomial coefficient in this range has a prime divisor in $2(N-1)$. This is a one-row assertion at both odd and even weights.

At any larger weight the same conclusion holds for an odd coefficient when $N$ is even and has at most six binary ones. The following more general result needs only a description of the parts. All assertions about larger weights here are subject to the stated hypotheses.

**Sparse-row theorem.** Suppose the three positive parts use at most six powers of an integer $B\ge2$ in total, counting repeated powers separately and allowing $B^0$. If $\gcd(B,N-1)=1$, then $K$ shares a prime with $N-1$. If $B$ is a prime $p$, the weaker condition $p^{v_p(N-1)}\le(N-1)^{3/5}$ suffices. In particular, for $p\mid N-1$ with base-$p$ digit sum $s_p(N)\le6$ and this size bound, **every** coefficient shares a prime with $N-1$. A coefficient avoiding $p$ has the required representation by Kummer's criterion. The same assertion follows if two different primes dividing $N-1$ both have digit sum at most six, since at most one of their full powers can exceed $(N-1)^{3/5}$.

## Two estimates for the part of the predecessor excluded from the gcd

Write $M=N-1$, $M_0$ for its largest divisor coprime to $B$, and $J=M_0/\gcd(M_0,K)$. Each full power of $J$ sees the parts as $0,0,1$ in some order. Indeed for $s^v\parallel M$, let $t=v_s(K)<v$. At each of the first $v$ layers the three residues sum to one or at least $1+s^j$ at modulus $s^j$. A failure to sum to one persists at all higher layers through $v$. Each such layer costs at least one in Legendre's formula. Thus the residues sum to one at $s^{v-t}$. Group the full powers of $J$ by the position with residue one. Call a group occupied when its product exceeds one.

- With at most one occupied group, $J$ divides two parts.
- With two occupied groups, $J<M^{2/5}$.
- With three occupied groups, $J<M^{12/35}$.

In particular if two rows of weight $N$ each satisfy the six-power hypothesis with a base coprime to $M$, and each has at least two occupied groups, their coefficients and the endpoint coefficient $N(N-1)$ share a divisor greater than $M^{1/5}$. Each coefficient contains a divisor of $M$ greater than $M^{3/5}$, and the least common multiple of those divisors divides $M$.

For the proof, a part consisting of one power of $B$ must be the exceptional position at every retained prime in $J$. Otherwise all three parts use exactly two powers. Write them $B^{i_j}(1+B^{d_j})$. In this case $M_0$ is odd. For odd $B$, the weight is a sum of six odd integers, and for even $B$ this follows by coprimality. If some $d_j=0$, that part again is coprime to $J$. It remains to treat three positive gaps with at least two occupied groups.

The [binomial gcd argument at a weight prime](six_units_weight_prime.md#retaining-divisors-of-the-endpoint) works for composite $B$ as well. Each occupied group divides the two factors $1+B^{d_j},1+B^{d_k}$ at its zero positions. A common odd prime forces $v_2(d_j)=v_2(d_k)$ by the order of $B$ at that prime. Two occupied pairs connect all three indices. Put $g=\gcd(d_1,d_2,d_3)$ and $t_{jk}=\gcd(d_j,d_k)/g$. Now the numerical gcd of the two plus factors is $B^{g t_{jk}}+1$. Bézout applied to the exponents gives $B^{g t_{jk}}\equiv-1$ modulo their gcd, because their normalized quotients are odd. The converse follows by factoring odd powers. The common factor $B^g+1$ divides $N$ and is coprime to $M$. The product of this group therefore divides
$$
 (B^{g t_{jk}}+1)/(B^g+1)<B^{g(t_{jk}-1)}.
$$
Every occupied index $t_{jk}$ is odd and at least three. The indices on different pairs are coprime. The two largest, $u<v$, share a vertex, so $g u v\le D=\max_j d_j$, and $B^D<M$.

For two occupied pairs, $(u+v-2)/(uv)\le2/5$ since $u\ge3,v\ge5$. For three, name them $a<u<v$. Then $u\ge5,\ a\le u-2,\ v\ge u+2$ and $(a+u+v-3)/(uv)\le12/35$. Multiplying the group bounds proves both estimates. The two integer inequalities are [checked in Lean](SparsePredecessorBins.lean).

If $K$ were coprime to $M$, **all three** exceptional positions would have to occur among the full powers of $M$. Otherwise a positive part would be divisible by $M$ despite being at most $N-2$. When $B$ is coprime to $M$, this gives $J=M$ with three groups, which is impossible. When $B=p$, removing the single full $p$-power leaves at least two groups, and $J=M/p^{v_p(M)}<M^{2/5}$ proves the prime-base assertion.

The verifier checks 39,160 rows, including 22,539 direct binomial products. As a positive control on the lowest congruences, the parts $178^3(1+178^{15}),178^2(1+178^{21}),178^3(1+178^{35})$ take the exceptional position successively at 29, 11 and 7. The actual coefficient valuations there are 29, 42 and 52.

## How the continuous range was checked

A row with odd coefficient consists of three disjoint binary submasks summing to $N$. If its coefficient is coprime to $M$, each full power of $M$ has a unique exceptional part. The products $A,B,C$ at these three positions are nontrivial, pairwise coprime, and $M=ABC$. For example the part at $A$ is the unique residue in $(0,M)$ that is one mod $A$ and zero mod $BC$.

The first program factors each $M$, constructs its idempotents by sums of the individual CRT residues, and tests them against the binary digits. Any possible row supplies six different admissible masks, the three parts and their three pairwise unions. A pair sum is less than $M$ because the omitted part exceeds one: a part equal to one would take every exceptional position. Test the remaining partitions by factorial floors.

The second program starts instead at the smallest of $A,B,C$, say $A$. Write $h=M/A>A^2$. For $h\equiv r\bmod A$ the first part is $a=h(r^{-1}\bmod A)$. Advance $a,N$ together by arithmetic progressions and test the binary submask condition. Split the bits of $N-a$ between $b<c$ directly. The conditions $M\mid ab$ and $M\mid b(b-1)$ force $a,b$ to be orthogonal idempotents modulo $M$, and leave $c\equiv1-a-b\bmod M$. Compute the exceptional groups of $b,c$ as $M/\gcd(M,b)$ and $M/\gcd(M,c)$ to check that $A$ is the smallest. Test the final rows by digit sums. This traversal needs neither a factorization table nor a CRT enumeration over all full powers.

| Inclusive weight bound | Odd rows passing the full-power congruences | Coefficient coprime to $N-1$ |
|---:|---:|---:|
| 20,000,000 | 128 | 0 |
| 1,000,000,000 | 548 | 0 |
| 2,000,000,000 | 684 | 0 |

The two enumerations agree at the two largest bounds. The [separate prime-pair calculation](odd_predecessor_prime_pairs.md) proves the same single-row divisibility at over six billion weights of the form $1+2^k p q$.

## Consequence for the equality in the part-gcd bound

Suppose one of the two companions to the endpoint row is $(g,gx,gy)$ with $x,y$ positive and coprime. Its coefficient is even, since at least one of $x,y$ is odd and its scaled part overlaps the lowest bit of $g$. Every prime dividing $N-1$ divides this coefficient too. Otherwise two parts would be zero modulo that prime. Since $g\mid N$, the two zero parts would have to be $gx,gy$, contradicting coprimality. Hence a second companion producing gcd one would have to be odd and coprime to all of $N-1$.

It follows, independently of checking other row shapes, that the equality $a=\gcd(b,c)$ in the [part-gcd theorem](part_gcd_gap.md) is excluded through two billion, at every even weight with at most six binary ones, and on the stated direct prime-pair range. At odd weights $10^9<N\le2\cdot10^9$ with precisely three predecessor primes the census also gives the full three-row assertion, since a counterexample requires one row avoiding all three, including 2. [Three-factor assignment](endpoint_factors.md)

~~~sh
mkdir -p .work
c++ -O3 -std=c++17 -pthread check_odd_predecessor_census.cpp -o .work/odd-census
c++ -O3 -std=c++17 -pthread check_odd_predecessor_parts.cpp -o .work/odd-parts
.work/odd-census 2000000000 8
.work/odd-parts 2000000000 8
python3 check_sparse_predecessor_bins.py
lean SparsePredecessorBins.lean
~~~
