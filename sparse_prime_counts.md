# Sparsity of exceptions at a fixed predecessor prime

**Theorem.** Fix a prime $p$ and an integer $s\ge7$. Among weights $N<p^{L+1}$ with $p\mid N-1$ and base-$p$ digit sum $s_p(N)=s$, at most
$$O_{p,s}(L^{s-4})$$
can fail Wasserman's three-row conjecture. For large $L$ there are $\Theta_{p,s}(L^{s-1})$ weights in this class, as is seen already by taking the $s-1$ positive digit positions distinct. Thus the proportion that could fail is $O_{p,s}(L^{-3})$. The estimate is ineffective. For $s=6$ the elementary [equal-gap theorem](equal_predecessor_gaps.md) gives $O(L^3)$ possible exceptions.

**Uniformity in the outer shift.** Normalize the $s-1$ positive positions by subtracting their minimum and use the other $s-2$ positions as coordinates. All patterns that can fail at *any* shift lie in a finite union of rational affine subspaces of codimension at least two. Away from finitely many inner spaces of codimension at least three, any one pattern of span $H>0$ can fail at only $O_{p,s}(1+L/H)$ shifts up to $L$. In particular, for seven units in any fixed prime base only $O_p(H^3)$ patterns of span at most $H$ can fail at any shift, out of $\Theta_p(H^5)$ possible patterns. For odd primes the [progression theorem](seven_progressions.md) gives an explicit exclusion through normalized span 20, with arbitrary dilation.

## Two disjoint sums of powers

Take disjoint sets $I,J$ of labelled coordinates of a projective torus with $t_0=1$ and $t_i=p^{v_i}$ for $1\le i\le n$. Each block has at least two coordinates, and at least one has three or more. Put $F_I=\sum_{i\in I}t_i$, and define $F_J$ likewise. For any $\delta>0$, away from a finite union of rational affine subspaces of codimension at least two in $\mathbb Q^n$,
$$\gcd_{(p)}(F_I,F_J)\le p^{\delta K},\qquad 0\le v_i\le K. \tag{1}$$
Here $\gcd_{(p)}$ removes every factor of $p$. It suffices to prove (1) for $K=\max v_i$ large, since a finite set of small exponent tuples fits in codimension two. The constants and exceptional spaces can depend on $p,\delta,I,J,n$.

First work on the full torus with $n+1$ variable coordinates. The common zero set $V$ of the two sums is irreducible of codimension two because it is a product of two nonempty torus opens in linear hyperplanes, and free torus factors. On a block of length $k\ge3$ no nontrivial character is constant on its zero set. Eliminate its last coordinate as minus the sum of the other $k-1$. In their rational function field these $k-1$ variables and their sum are pairwise nonassociate irreducibles, so all their powers in a constant monomial must be zero. On a block of length two the constant characters are exactly the powers of the ratio of the two coordinates, and that ratio equals $-1$. Free coordinates have exponent zero in a constant character. Since the blocks vary independently, the only possible nontrivial primitive characters constant on $V$ are the ratio in the length-two block and its inverse, both with value $-1$.

It follows that the two sums stay coprime in the Laurent polynomial ring on every connected codimension-one torus translate containing a positive point. In fact a common divisor there would give a common zero locus of dimension at least $n-1$, equal to $\dim V$. Irreducibility would force $V$ to lie in the translate. A connected codimension-one translate is cut out by a primitive character with a constant value. The only candidates have value $-1$ on $V$, whereas every character is positive at a positive point. The same coprimality holds after fixing $t_0=1$. Indeed the full torus is the product of this normalized torus with a scalar coordinate. Both sums are homogeneous of degree one, so common factors on a normalized translate pull back to common factors on its full inverse image. Coprimality on the ambient normalized torus follows from the codimension of $V$ in the same way.

Apply the outside-$S$ gcd theorem with places $S=\{\infty,p\}$ to the two normalized sums. For rational integer values the outside-$S$ logarithmic gcd is the logarithm of $\gcd_{(p)}$. The theorem gives (1) away from a finite union of proper torus translates. Split into connected components and discard components containing no exponent point. On the points considered, independent character equations become independent rational affine equations on $v$. This already handles every component of codimension at least two. [Levin–Wang, Theorem 1.13(a), stating the Corvaja–Zannier and Levin gcd theorem](https://doi.org/10.1515/crelle-2019-0033)

For each codimension-one component meeting the points, take an integer solution of its primitive affine equation and a lattice basis of the homogeneous solutions. This writes its points as $v=v_0+M w$, with $\|w\|_\infty=O(1+K)$, and parametrizes the connected torus translate by the corresponding Laurent monomials. The two restricted sums are coprime by the preceding argument. Clear coordinate monomials and powers of $p$ to obtain ordinary coprime polynomials, leaving all valuations away from $p$ unchanged. A Laurent unit contributes only a fixed factor to the outside-$S$ gcd. Otherwise apply the same theorem in the coordinates $p^{w_j}$, with an exponent small enough for the fixed height comparison. Its exceptional translates now have codimension at least two in the original torus. Taking the finite union proves (1).

The same proof permits nonzero rational coefficients in a block of length at least three. If one block has length two, its two coefficients must have the same sign, so that their ratio with a minus sign is negative. Clearing fixed coefficient denominators changes the logarithmic bound by at most a constant, absorbed by first taking a smaller exponent.

The codimension in this gcd statement can be attained. With $t=p^k$ the sums $p^j(1+t^3)$ and $1+t^2+t^4$ have outside-$p$ gcd exactly $t^2-t+1$. Their projective exponent vectors $(j,j+3k,2k,4k)$ form a two-dimensional set in four dimensions. The identity follows by removing the common factor and using $\gcd(t+1,t^2+t+1)=1$.

## Three disjoint sums

For three disjoint blocks, each of size at least two, suppose at least one has size at least three. For positive unit coefficients the outside-$p$ gcd of the three sums obeys (1) off a finite union of subspaces of **codimension at least three**. More generally, nonzero rational coefficients are allowed provided that if two blocks have size two, at least one of these two has coefficients of the same sign.

In the full torus the common zero set is now irreducible of codimension three, by the same product description. On a codimension-one translate the restricted sums cannot have a common nonunit factor, since such a factor would have a zero locus of codimension two in the ambient torus. On a connected codimension-two translate a common factor would force the whole joint zero set to lie in the translate, by equality of dimensions. This requires two independent constant characters. The only possibility is that exactly two blocks have length two. Their two ratios generate a saturated rank-two lattice of characters. The rank-two lattice defining a connected translate is also saturated, so if it is contained in this lattice it must equal it. In particular *each ratio separately* would be fixed on the translate. At least one ratio has negative value on the joint zero set, excluding every translate containing a positive point.

To apply the gcd theorem whenever three Laurent polynomials have no common nonunit factor, choose one nonzero member and an integer linear combination of the other two avoiding every irreducible factor of the first. Each such factor excludes only a proper linear subspace of coefficient choices. The outside-$p$ gcd of the three divides that of the chosen pair at integer evaluations, and the same valuation inequality holds at rational evaluations with powers of $p$ in their denominators. A Laurent unit gives a bounded outside-$p$ gcd directly. Apply the gcd theorem first on the whole torus, then on the relevant components of codimension one, and finally on any resulting components of codimension two, using the same lattice parametrization and height comparison as above. The remaining components have codimension at least three.

## The factors assigned to two rows

Write the positive positions as $a+h_i$, $1\le i\le s-1$, with $a\ge1$ and $0=h_1\le\cdots\le h_{s-1}=H$. Identical positions are separate labels. Put $W=p^a$ and $m=\sum p^{h_i}$, so $N=1+Wm$, $p\nmid m$ and $p^H\le m\le(s-1)p^H$. In a counterexample select a nonendpoint row avoiding $p$ and assign $p$ to it. The [endpoint theorem](all_weights_endpoint_reduction.md) and [prime-power assignment](endpoint_reduction.md) give its parts as
$$(Wx,Wy,1+Wz),\qquad x+y+z=m,$$
where Kummer's criterion partitions the labelled powers among $x,y,z$. Its assigned weight factor $U\ge3$ divides all three parts. Thus $x,y$ each use at least two units and $z$ at least one. The complementary weight factor is $V=N/U$. If $S$ is the product of full prime powers of $m$ assigned to the first row, each such power divides $x$ or $y$, since two actual parts vanish modulo it. Each remaining power of $m$ divides two normalized parts of the companion row after division by $V$. Those three normalized parts sum to $U$. Consequently
$$S\mid\gcd(m,x)\gcd(m,y),\qquad 27m^2\le U^3 S^2. \tag{2}$$
Both $S,U$ are coprime to $p$.

Choose $\delta=1/10$. Split $S$ into three factors according to the two zero positions, assigning each full power to its pair of zero positions. Write $S_{xy}$ for the factor on the two unshifted parts, and $S_x,S_y$ for the factors on the shifted part and respectively $Wx,Wy$. Then
$$S_{xy}\mid D=\gcd(x,y,z),\qquad
S_x\mid\gcd(x,m),\qquad S_y\mid\gcd(y,m),\qquad
U S_x S_y\mid 1+p^a z. \tag{3}$$
For the first assertion use $m=x+y+z$ and $S_{xy}\mid m$. The three factors on the last left side are pairwise coprime.

If one of $x,y,z$ is a singleton then $D=1$. If all use at least two units and one uses at least three, the three-block version of (1) gives $D\le m^\delta$ off a finite union of codimension-three inner spaces. Here $D$ is coprime to $p$ since it divides $m$. A rational affine space of dimension $k$ contains $O((1+L)^k)$ integer points in a box of side $L$, by projecting injectively to some $k$ coordinates. Thus the exceptional spaces in these $s-2$ inner coordinates contribute only $O_{p,s}(L^{s-5})$ patterns times at most $L$ shifts.

There is one remaining possibility for $D$. If $x,y,z$ are all pairs, a zero gap gives $D\le2$. For three positive gaps, a $D>2$ divides $p^g+1$ where $g$ is their gcd and all three quotients by $g$ are odd. Indeed Bézout gives $(p^g)^2=1\bmod D$, and each quotient must then be odd to give $-1$ at the three gaps. If all gaps coincide, $U\mid p^g+1\mid m$ is impossible. Otherwise $g\le H/3$, whence $D\le p^{H/3}+1$.

It follows in every remaining case, for sufficiently large $H$, that at least one of $U,S_x,S_y$ exceeds $m^\delta$. Otherwise (2) and (3) would give $27m^2\le m^{7\delta}D^2$, impossible even with $D\le p^{H/3}+1$. Bounded $H$ contributes $O_{p,s}(L)$ weights.

## How many outer shifts can supply a large divisor?

We use the following elementary estimate. Fix $c,\delta>0$, $K\le p^{cH}$ and a positive integer $z$. There are only $O_{c,\delta}(1+L/H)$ integers $1\le a\le L$ with
$$\gcd(K,1+p^a z)\ge p^{\delta H}. \tag{4}$$
Indeed choose a fixed integer $k\ge2$ with $k\delta>c$, then $\eta>0$ with $k\delta-\binom{k}{2}\eta>c$. If $k$ such shifts lie in an interval of length $\eta H$, let $D_j$ be their individual gcds. Each pairwise gcd divides $p^{|a_i-a_j|}-1$ by subtraction, since it is coprime to $pz$. Prime valuations give
$$K\ \ge\ \operatorname{lcm}_j D_j\ \ge\
\frac{\prod_j D_j}{\prod_{i<j}\gcd(D_i,D_j)}
\ \ge\ p^{(k\delta-\binom{k}{2}\eta)H}>K.$$
Partitioning the interval of shifts proves (4).

Suppose the large factor is $S_x$ or $S_y$. For $\gcd(m,x)=\gcd(x,y+z)$ the two disjoint blocks have lengths at least two and three. The same holds for $\gcd(m,y)$. Thus (1), with exponent $\delta/2$, places the pattern in a finite union of inner spaces of codimension at least two. If $U$ is large and either $x$ or $y$ uses at least three units, use $U\mid\gcd(x,y)$ to reach the same conclusion. Moreover, whenever any one of $U,S_x,S_y$ exceeds $m^\delta$, the integer $K=m\gcd(x,y)\le m^2\le p^{3H}$ for large $H$ is a multiple of that factor. By (3) the shifts satisfy (4). This proves the asserted bound on the number of shifts per pattern outside the codimension-three inner spaces, including the case of two unshifted pairs. At $s=7$ every allocation either has a singleton in $z$ or has three inner pairs. The argument therefore bounds the number of shifts for every seven-unit pattern, without any exceptional inner spaces.

## Reducing the shift when two unshifted pairs have large gcd

It remains to count $U>m^\delta$ when both $x,y$ are pairs. Write $x=p^i(1+p^d)$ with $d>0$ since $U\ge3$ and $p\nmid U$. Then $p^d=-1\bmod U$ and $U\le p^d+1$, so $d\ge\delta H/2$ for sufficiently large $H$. Take $0\le a'<2d$ congruent to $a$ modulo $2d$. The three sums
$$p^{a'}x,\qquad p^{a'}y,\qquad 1+p^{a'}z$$
still vanish modulo $U$. They consist of disjoint blocks of lengths $2,2,s-4$ on all $s$ projective coordinates, the last block having at least three terms. Their $s-1$ variable exponents are $a'+h_i$, all lying between zero and $3H$. By the three-block lemma with exponent $\delta/4$, a gcd greater than $m^\delta$ puts this reduced vector in a finite union of spaces of dimension at most $s-4$.

The minimum of the reduced vector is its first coordinate $a'$, since the positions were sorted. Subtracting it projects these spaces linearly into the inner position space of dimension $s-2$. Their images have dimension at most $s-4$. Together with the previously obtained inner spaces and finitely many bounded patterns this proves the assertion on patterns, uniformly in the outer shift.

For $R\le H<2R$ outside the codimension-three inner spaces there are consequently $O_{p,s}(R^{s-4})$ possible patterns, each admitting at most $O_{p,s}(1+L/R)$ shifts. Summing over dyadic $R\le L$ gives $O_{p,s}(L^{s-4})$, since $s\ge7$. Including the codimension-three contribution already counted and the finitely many label assignments completes the proof.

## When the sparse prime divides the weight

A related statement applies at $p\mid N$. Fix $p$, $s\ge7$ and $0<\rho<2$. Write $a=v_p(N)$ and $H$ for the span of the $s$ labelled base-$p$ digit positions, so their minimum is $a$. Among $N<p^{L+1}$ with $s_p(N)=s$ and $a\le(2-\rho)H$, at most $O_{p,s,\rho}(L^{s-3})$ weights can fail. More precisely their $s$ exponents lie in a finite union of rational affine subspaces of codimension at least three.

Assign $p$ to a row avoiding it. Write $W=p^a$, $N=Wm$ and the row as $(Wx,Wy,Wz)$, where $m=x+y+z$ and the parts partition the $s$ digit units. The weight factor assigned to the row is $Wu$, with $p\nmid u$ and $u\mid\gcd(x,y,z)$. For each full power of $N-1$ assigned to this row, choose its two zero positions and denote the products of these powers by $S_{12},S_{13},S_{23}$. Thus, for instance,
$$S_{12}\mid\gcd(Wx,Wy,Wz-1).$$
For their product $S$ the companion product bound reads $27(N-1)^2\le (Wu)^3S^2$. All four factors $u,S_{12},S_{13},S_{23}$ are coprime to $p$.

Use the original exponent coordinates and adjoin the constant projective coordinate 1. For $u$ take the three sums $Wx,Wy,Wz$. For $S_{12}$ take the first two unchanged and the third with the additional term $-1$, and permute for the other choices. In each instance the three blocks are disjoint. A singleton among the sums without a constant term gives outside-$p$ gcd 1 immediately. Otherwise all blocks have size at least two and at least one has size at least three, since there are $s$ terms for $u$ and $s+1$ for each $S_{ij}$. If a block of size two uses the negative sign, the other block of size two, if present, has both signs positive. The three-block lemma therefore applies to each gcd.

The maximum exponent is $a+H\le3H$. Taking $\delta=\rho/54$, apart from the stated codimension-three spaces each of the four factors is at most $p^{3\delta H}$ for large $H$. It would follow that
$$27(N-1)^2/W^2\le Wp^{27\delta H}\le p^{(2-\rho/2)H}.$$
But $N-1\ge Wm/2$ and $m\ge p^H$. This proves the assertion because bounded $H$ also bounds $a$ in this range.

In particular, among *all* multiples of a fixed $p$ of fixed digit sum $s\ge7$, the upper limiting proportion of exceptions, measured through $p^{L+1}-1$, is at most $3^{1-s}$. To see this, repeated exponents account for $O_s(L^{s-1})$ of the $\Theta_s(L^s)$ weights. Conditional on maximum exponent $j$ and distinct positions, the other $s-1$ positions form a subset of $\{1,\ldots,j-1\}$. If $a>(2-\rho)H$, all lie above $j(2-\rho)/(3-\rho)$. The proportion tends to $(3-\rho)^{1-s}$ as $j$ tends to infinity, so the same holds after summing over $j\le L$. On the complementary range the exceptions have zero limiting proportion by the proved bound. Let $\rho$ decrease to zero.

A numerical check of the underlying integer identities, the reduced period and the large-divisor spacing estimate is provided in the [arithmetic verifier](check_sparse_prime_counts.py). It checks 40,500 period reductions, including 6,750 reductions to exponent zero where only the outside-$p$ gcd is invariant. It also tests 14,336 shifts with deliberately constructed large divisors and 1,650 instances of the cyclotomic gcd example. The outside-$S$ gcd theorem supplies the affine exceptional spaces nonconstructively.
