# Counting pairs of primitive directions across all weights

Let $X\ge3$. Choose two primitive positive integer triples whose respective sums are at most $X$, then scale them to any common weight $N$. **Only $O_\epsilon(X^{5/2+\epsilon})$ pairs can pass the lowest-layer conditions at even one common weight.** The number of direction pairs before these tests is of order $X^6$, by the positive density of primitive triples. Every excluded pair has a common coefficient prime dividing $N(N-1)$ at every scale. This count uses no high-digit carry tests.

More precisely, if $g$ is the gcd of all six actual entries together, the number of pairs that can pass at some $g\ge Z$ is
$$
              O_\epsilon\!\left(X^\epsilon
                     (X^2+X^{5/2}/Z^2)\right),\qquad Z\ge1.
$$
In particular, for scales at least $X^{1/4}$ the count is $O_\epsilon(X^{2+\epsilon})$. Each surviving pair has at most one possible weight by the [complementary-divisor formula](primitive_row_scales.md).

### Counting by an assigned divisor

Let $x,y$ have respective sums $u,v$. In a surviving pair the sums are coprime, $N=guv$, and the full powers of $N-1$ can be assigned to the two directions. Call their products $S,T$. At each power in $S$ two entries of $x$ vanish and the third is $u$ modulo the power, and similarly for $y,T,v$. In particular $ST=guv-1$. One row has assigned product at least its sum, since otherwise $ST\le(u-1)(v-1)$. Name it $x$, so $S\ge u$.

Work in dyadic ranges for $u,v,S,T$, and denote their sizes by $U,V,D,E$. Split into the cases below using the actual integers before applying the dyadic estimates. The comparison symbols absorb absolute dyadic constants. The product of the three entries of $x$ is divisible by $S^2$, hence $S\le(u/3)^{3/2}$, and likewise $T\le(v/3)^{3/2}$. In particular $ST+1\ll X^3$.

Group the powers in $S$ by the exceptional position, with group products $d_1,d_2,d_3$. Then $d_1 d_2 d_3=S$ and $S/d_i$ divides $x_i$. For given groups the number of such positive triples with sum in the dyadic range is at most $O(U^3/S^2)$, by counting multiples independently. There are $O_\epsilon(X^\epsilon)$ groupings by the divisor bound. Summing over $S$ in its dyadic interval therefore gives $O_\epsilon(X^\epsilon U^3/D)$ possibilities for $(x,S)$.

Now $ST\equiv-1\pmod u$ and $\gcd(S,u)=1$. There are $O(1+E/U)$ choices for $T$. Then $v\mid (ST+1)/u$ has $O_\epsilon(X^\epsilon)$ possibilities. At every power of $T$, the exceptional position of $y$ fixes the three residues modulo that power. After the $O_\epsilon(X^\epsilon)$ possible assignments, the first two coordinates are each fixed modulo $T$ and the sum determines the third. The count in the box, with the divisor exponents absorbed into an arbitrarily small $\epsilon$, is consequently
$$
       \ll_\epsilon X^\epsilon\frac{U^3}{D}
                       (1+E/U)(1+V/E)^2. \tag{1}
$$

If $T\ge v$ as well, order the two sides so $S\ge T$. Both assigned products reach their respective sums. Then (1) is $O_\epsilon(X^\epsilon U^2)$ since $D\gg U,E$.

Suppose $T<v$. Use $DE\gg ZUV$. For $T\ge u$, (1) is at most $O_\epsilon(X^\epsilon U^2V^2/(DE))$, which is $O_\epsilon(X^\epsilon UV/Z)$. For $T<u$ it is at most $O_\epsilon(X^\epsilon U^2V/(ZE))$. Since $D\ll U^{3/2}$, we have $E\gg ZV/\sqrt U$, yielding $O_\epsilon(X^\epsilon U^{5/2}/Z^2)$. Summing the logarithmically many boxes and readjusting $\epsilon$ proves the result.

### A finite census independent of a weight cutoff

The [enumerator](check_direction_pair_count.py) generates the grouped divisors of the first row directly. Each $d_i>1$ because a positive entry cannot be a multiple of $S\ge u$. They are pairwise coprime, so up to permutation take $2\le d<e<f$. The entries are positive multiples of $ef,df,de$ with sum at most $\min(X,def)$. The program retains primitive rows, enumerates $T$ in its class modulo $u$, then obtains all companion rows in the residue classes from their assigned powers. It validates each surviving pair by the independent scale formula and checks the actual coefficients by Legendre valuations without constructing huge factorials.

As a cross-check, through primitive sum 200 it also enumerates the primitive triples by their entries and groups them by their sum and the gcd of their pair products. For each possible pair of profiles it tries the divisors of the lcm of these gcds directly. The two methods recover the same eight surviving pairs.

Through primitive sum **6,000**, with no cutoff on the common weight, the run generated 1,942,955 first rows with assigned modulus. Exactly **551,062** unordered primitive direction pairs passed the lowest layers, at 127,495 distinct weights. Each pair had a common endpoint prime under the full valuation tests. The largest surviving weight was **1,548,955,976**. The smaller direct cross-check covers 17,327,493,561 unordered pairs through sum 200 by grouping profiles. Run `python3 check_direction_pair_count.py --limit 6000` with Python 3.9 or later.

Among the 551,062 surviving pairs the smallest endpoint gcd was **42**, attained at $N=1550578$ by directions $(16,1571,1571)$ and $(21,141,329)$. Every surviving pair shared at least two endpoint primes. The gcd **44** at $N=8055653$, with directions $(44,892,2453)$ and $(735,821,821)$, attains that minimum count. Adding `--exact-records` with python-flint installed recomputes both gcds with exact binomial products. The unique surviving pair at the largest weight has directions $(791,2100,2825)$ and $(1073,2117,2701)$.
