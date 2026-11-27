# Powers of two near prescribed trinomial proportions

Write $K(x_1,x_2,x_3)=N!/(x_1!x_2!x_3!)$ for positive integer parts with sum $N$, and call $(x_1/N,x_2/N,x_3/N)$ their proportions.

**For area-almost every point of the open triangle of proportions and every $J\ge1$, there is a neighborhood in which every trinomial coefficient at every weight is divisible by $2^J$.** In particular, for almost every choice of three limiting proportions, if three rows of the same weight tend to their respective limits, their common divisor is divisible by arbitrarily high powers of two. The neighborhood can depend on the point and on $J$.

There is a conditional result along individual line segments, which themselves have zero area. Fix odd positive integers $h<d$. For length-almost every $\beta\in(0,1-h/d)$:

* Near $(h/d,\beta,1-h/d-\beta)$ every coefficient is even.
* For each $J\ge1$ there is a neighborhood of that point in which $2^J$ divides every coefficient whose first proportion is **at most** $h/d$.

The second assertion is one-sided for a reason. At $N=5\cdot2^m-1$ take $x_1=3\cdot2^m$ and split $2^{m+1}-1$ into two arbitrary positive parts. Those two parts add without a binary carry, since they are complementary masks. Adding their sum to $x_1$ gives exactly two carries, at positions $m$ and $m+1$. Thus every point of the open segment at first proportion $3/5$ is approached from above with valuation exactly two.

## Radial projections of digit sets

Here is the geometric argument used below, in ambient dimension $l+1$ for $l=1,2$. Suppose a finite collection of compact sets $E_i$ has the following two properties.

1. Each $E_i$ at level $n$ is covered by $O(2^{ln})$ translated copies of $2^{-n} E_j$, with no rotations and a constant independent of $n$.
2. All coordinate $l$-plane projections of every $E_i$ have zero $l$-dimensional measure.

Then their radial projections from **any** point have zero $l$-dimensional measure on the unit sphere, omitting the center. We give the argument also for the finite graph of digit states used here. Simon and Solomyak prove the planar self-similar visibility theorem, including systems with rotations. [Simon and Solomyak, Theorem 1.1](https://arxiv.org/pdf/math/0510342)

Property 1 gives finite $\mathcal H^l$ measure by the defining covers. Property 2 implies pure $l$-unrectifiability. On a countably rectifiable $l$-set its approximate tangent exists almost everywhere. Some coordinate projection has positive tangential Jacobian at every tangent plane, since some maximal minor is nonzero. The area formula for these finitely many projections rules out a positive-measure rectifiable part of $E_i$. By the Besicovitch–Federer projection theorem, the orthogonal projection of each $E_i$ to $u^\perp$ has zero $l$-measure for almost every direction $u$. The identification of directions up to sign with codimension-one planes preserves the rotation-invariant measure classes. [Federer, *Geometric Measure Theory*, area formula and Theorem 3.3.13](https://link.springer.com/book/10.1007/978-3-642-62010-2)

To pass to radial projection, translate the center to zero, restrict to an annulus away from zero and work in a smooth coordinate chart on the sphere, on a compact subpatch. For the charted direction map $F$ the derivative at $x$ factors through orthogonal projection onto $(x/|x|)^\perp$ and a uniformly bounded linear map. Fix a compact subset $C$ of good directions in the subpatch. The $l$-volumes of the $\eta$-neighborhoods of the projections of all $E_j$ tend uniformly to zero for $u\in C$ as $\eta\downarrow0$. To see this, each compact null projection has a neighborhood of arbitrarily small volume, and nearby projection planes can be identified by rotations close to the identity. Compactness supplies a finite cover of $C$.

For a level-$n$ copy meeting the preimage of $C$, linearize $F$ at one such point. On that copy the remainder is $O(2^{-2n})$ by a uniform bound on second derivatives, using a slightly larger annulus and chart. Its image in the patch thus has volume $2^{-ln}o(1)$ uniformly, by the preceding projected-neighborhood estimate. Equivalently one can apply the estimate to the linear image in chart coordinates, where the map on the projected plane and its inverse have bounded norm on the compact sets in question. Summing property 1 gives zero volume over $C$. Inner regularity, the null set of exceptional directions, and a countable exhaustion by patches and annuli give the result. The same conclusion holds for any countable union of translated and dilated copies of these sets.

## All three proportions free

Let $E$ consist of three real numbers in $[0,1]$ given by infinite binary expansions whose joint digit in each column is one of $000,100,010,001$. It has a cover by $4^n$ translated copies at scale $2^{-n}$. Each coordinate-plane projection is covered by $3^n$ squares of side $2^{-n}$, so has area zero. The preceding argument applies with $l=2$. Passing from the sphere to the sum-one triangle is smooth with a bounded derivative in the nonnegative octant.

For a row with $v_2(K)<J$, at most $J-1$ binary columns have two or more ones among the parts. Each such column sends a carry to the next place, and Kummer's formula counts all outgoing carry amounts with multiplicity. Scale the parts by the same binary power to make their sum lie in $[1/2,1]$. Pad their digit words with zeros. The compact space of infinite words with at most $J-1$ such columns includes all the normalized rows. Every one of its words eventually has a tail in $E$. Thus its geometric image belongs to a countable union of translated and scaled copies of $E$. On restricting the sum to $[1/2,1]$ its directions form a **compact null set** by the radial projection argument. Outside the union over $J$, take the complement of each compact set as the required open neighborhood.

## An odd rational first proportion

Consider normalized rows approaching first proportion $h/d$ with bounded binary valuation. Pass to a convergent subsequence of the padded digit words of the total and the three parts. Also retain at each boundary the incoming carry from the lower digits. The carry amounts are at most two and their sum is bounded. The limiting stream therefore has zero carries at and beyond some depth. At such depths the prefix $n$ of the total equals the sum of the three part prefixes. If $a$ denotes the first part prefix, put $z=da-hn$.

The real coordinates in the limit satisfy $d x_1=h(x_1+x_2+x_3)$. The three fractional tails after a sufficiently deep prefix have sum at most one and give
$$
                 h-d\le z\le h .
$$
Indeed $z=-(d-h)u_1+h(u_2+u_3)$ for tails $u_i\in[0,1]$ with $\sum u_i\le1$. Under the one-sided assumption $dx_1^{\rm int}\le hN$ on the original integer rows, $da<h(n+1)$ at **every** prefix of each row, since $N<(n+1)M$ and $aM\le x_1^{\rm int}$ at scale $M$. The strict integer bound passes to the limit. If the coefficients are odd, no one-sided assumption is needed. All prefixes have zero incoming carry, the same closed bounds apply from the start, and the state $h$ cannot first be reached. Indeed the transitions are
$$
     z\mapsto 2z,\qquad 2z+d-h,\qquad 2z-h,\qquad 2z-h .
$$
As $d,h$ are odd, a transition into the odd state $h$ must be of the last type and start at $h$. The initial state is zero.

Use the $d$ states $h-d,\ldots,h-1$. Keep transitions remaining in that range. They correspond respectively to the new bit being absent or belonging to the first, second or third part. Label each edge by the binary digits of the total and the second part. Each target state has at most two incoming edges, counting multiplicities. For an even target only the first two possibilities occur. For an odd target only the last two occur, with the same predecessor. Thus there are at most $d\,2^n$ paths of length $n$, and the compact sets of two real coordinates generated from each state satisfy property 1 with $l=1$.

Both coordinate projections have length zero. Choose $L$ with $2^L>d$. A run of $L$ zero total digits must start at state zero, since it multiplies the starting state by $2^L$. A single occupied digit from zero cannot return to zero. Thus $0^L1 0^L$ is forbidden in the total. In the second part a run of $L$ ones sends $z$ to $h+2^L(z-h)<h-d$ and is forbidden. Avoiding any fixed binary word gives a null set, by partition into disjoint blocks of its length. This proves property 2.

After the finitely many exceptional columns each limiting pair (total, second part) consequently lies in a translated and scaled graph set, with the earlier digits providing the translation. The ratio of the second coordinate to the first has zero length, since this ratio is a smooth coordinate on the radial directions in the normalized range. The relevant limiting proportions form a compact set: take the closure of the original directions of bounded valuation, with the one-sided restriction when applicable, and intersect with the line. The subsequence argument covers every point of that intersection and also includes individual rows by zero padding. This proves the line-segment assertions.

## Whole progressions of weights at cofactor five

For an odd prime power $Q$ consider $N=1+5\cdot2^kQ$. Call a residue class of $k$ certified when the full conjecture holds for **every positive** $k$ in that class. The following assertions refer to the entire period $\operatorname{ord}_Q(2)$.

* For each fixed odd prime $\ell$, as $e\to\infty$ the fraction of certified classes at $Q=\ell^e$ tends to one. The same holds within any fixed arithmetic progression of $k$, taking its induced period.
* For every $\epsilon,\delta>0$, all but $O_{\epsilon,\delta}(X^\delta)$ odd primes $Q\le X$ have at least a fraction $1-\epsilon$ of their classes certified. This also holds within any fixed progression, with its step allowed in the implied constant.

Here one can select classes without knowing any counterexample. The [cofactor-five reduction](five_leading_digits.md) requires the odd row $(3PQ,cQ,(2P-c)Q+1)$, where $P=2^k$ and $c$ is the least positive inverse of $Q$ modulo $2P$. Put $j=(cQ-1)/(2P)$, so $0\le j<Q$ and $j\equiv-2^{-k-1}\bmod Q$. The first share approaches $3/5$ from below uniformly in $k$ as $Q\to\infty$. The second share differs from $2j/(5Q)$ by $O(1/Q)$ uniformly. The null compact set on the $3/5$ segment implies that given any $\eta>0$, there is a fixed finite union $U$ of intervals of total length less than $\eta$ such that every failing pair $(k,Q)$ with $Q$ sufficiently large satisfies $j/Q\in U$. Take a small open cover of the compact null set after scaling by $5/2$ and use sequential compactness. Such covers can also be found by finite binary digit search. The intervals (1) in the [leading-digit calculation](five_leading_digits.md#projecting-the-leading-cells-onto-the-fraction-r) give nested covers at growing binary depths. Every point surviving all depths has an infinite path by finite branching, and the prefix interval widths tend uniformly to zero. It thus belongs to the null ratio set. Finite covers of arbitrarily small total length result. The $O(1/Q)$ discrepancy can be included by slightly enlarging the chosen intervals.

For $k=k_0+s n$ with fixed $s\ge1$, let $g=2^s$, $L=\operatorname{ord}_{\ell}(g)$ and $v=v_{\ell}(g^L-1)$. For $e>v$, on each fixed class of $n$ modulo $L$ the values of $j$ fill all $\ell^{e-v}$ lifts of one unit class modulo $\ell^v$. Indeed raising $1+u\ell^i$ with $i\ge1$ and $\ell\nmid u$ to the odd prime power exponent $\ell$ increases its depth by exactly one, by binomial expansion. The cyclic subgroup generated by $g^L$ therefore has precisely the size of that group of principal units. The normalized $j$ form a translate of an equally spaced grid. If $U$ has $m$ intervals, the fraction it captures is at most $|U|+O(m/\ell^{e-v})$. This proves the first assertion. It is uniform across bases whenever the displayed number of lifts tends to infinity.

For $Q$ prime, the $j$ constitute a coset of the subgroup generated by $2^s$ modulo $Q$. Fix $0<\alpha<1$ with $2\alpha<\delta$. A subgroup of size at least $Q^\alpha$ has all normalized nontrivial additive Fourier coefficients $O(Q^{-c(\alpha)})$ by the Bourgain–Glibichuk–Konyagin subgroup estimate, uniformly for multiplicative cosets. Thus the fraction in the fixed finite union $U$ is at most $|U|+o(1)$ uniformly in such primes and cosets, by approximating its indicator from above by trigonometric polynomials. All primes through $X$ with smaller order divide $\prod_{i\le X^\alpha}(2^{si}-1)$, whose logarithm is $O_s(X^{2\alpha})$. This bounds even the number of its distinct prime factors by $O_s(X^{2\alpha})$. Taking $\eta<\epsilon$ proves the second assertion. [Bourgain, Glibichuk and Konyagin, subgroup exponential-sum theorem](https://doi.org/10.1112/S0024610706022721)
