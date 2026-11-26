# Effective exceptional exponent spaces, uniform in the prime

Fix an integer $s\ge7$. Consider weights with $p\mid N-1$ and $s_p(N)=s$, where $p$ is any prime and the digit sum counts units with multiplicity. Write the positive digit positions as $a+h_i$ with $a\ge1$ and $0=h_1\le\cdots\le h_{s-1}=H$.

**Inner positions.** There is an effectively computable finite union of rational affine subspaces in the $s-2$ coordinates $(h_2,\ldots,h_{s-1})$, each of codimension at least two, that contains every pattern which can fail Wasserman's conjecture for **any prime and any outer shift**. The same conclusion holds for patterns which can fail at some positive integer dilation. In particular at most $O_s((1+L)^{s-4})$ inner patterns of span at most $L$ can fail at any prime, shift or dilation.

**Full positions.** For each prime the number of possible exceptional weights below $p^{L+1}$ in this class is $O_s((1+L)^{s-4})$, with an effective constant independent of $p$. If $a\le C H$ for fixed positive rational $C$, the full $s-1$ exponent tuples of all possible exceptions, **combined across primes**, lie in an effectively computable finite union of rational affine spaces of codimension at least three.

For example, take $n=s-2$ and $h_1=0$. Use the other $n$ positions
$$b\big((i+1)d^n+d^i\big),\qquad 0\le i<n.$$
The conjecture holds for these positions simultaneously for **all** primes, shifts $a\ge1$ and dilations $b\ge1$ whenever $d\ge 2n(1+\lceil\log_2(n+1)\rceil)+1$. At **seven units** ($s=7$) it holds at **every integer $d\ge2$**. The positions become nearly equally spaced as $d$ grows.

The [previous sparse-position estimate](sparse_prime_counts.md) used an outside-$S$ gcd theorem for a fixed $p$. The effective argument here uses integer determinants and simultaneous rational approximation of the exponents.

## An elementary common-base gcd bound

Let $F,G\in\mathbb Z[X_1,\ldots,X_n]$ be nonzero and coprime in the Laurent polynomial ring over $\mathbb C$. Allow their coefficients to depend on the prime $p$. Suppose their total degrees are at most $D\ge1$ and their coefficient $1$-norms at most $p^B$, with fixed nonnegative integer $B$. Let $e\in\mathbb Z^n$, $K=\|e\|_\infty$, and write $p^e=(p^{e_1},\ldots,p^{e_n})$. Divisibility by a positive integer $W$ coprime to $p$ is taken in $\mathbb Z[1/p]$.

For every rational $\epsilon>0$, if $W$ divides both $F(p^e)$ and $G(p^e)$ and $W>p^{\epsilon K}$, then $e$ lies in an effectively given finite union of proper rational affine hyperplanes. The list depends only on $n,D,B,\epsilon$. In particular the same list works even if the polynomial pair varies within the stated hypotheses.

Choose $Q=\max(2,\lceil16D^2/\epsilon\rceil)$ and put $q_0=Q^{n-1}$. For $K>0$ there are $1\le q\le q_0$ and $v\in\mathbb Z^n$ such that $\|v\|_\infty=q$ and $\|q e/K-v\|_\infty\le1/Q$. One coordinate of $e/K$ is $\pm1$. Apply the pigeonhole principle to the fractional parts of the first $q_0+1$ multiples of the other coordinates in $Q^{n-1}$ boxes. The integer approximation in each coordinate has absolute value at most $q$.

Give each variable $X_i$ weight $v_i$. For $F=\sum f_\nu X^\nu$, let $l_F$ be its least monomial weight, and form
$$F_v(T)=\sum_\nu f_\nu X^\nu T^{\nu\cdot v-l_F}.$$
Do the same for $G$. If either polynomial has degree zero in $T$, use its coefficient as $R(X)$. Otherwise let $R=\operatorname{Res}_T(F_v,G_v)$. This is nonzero. In the Laurent polynomial ring with $T$ adjoined, substituting $X_iT^{v_i}$ for $X_i$ is invertible and preserves coprimality. Localizing the coefficient variables to their rational function field also preserves absence of common factors, and neither shifted polynomial has a factor $T$. At $T=1$ the resultant identity gives $W\mid R(p^e)$.

All monomials of $R$ have the same $v$-weight. Indeed under $X_i\mapsto X_i t^{v_i}$, the polynomial $F_v(T)$ becomes $t^{l_F}F_v(tT)$, and similarly for $G_v$. Scaling the two polynomials and the argument multiplies their resultant by one power of $t$. This also holds for the degree-zero choice by its definition.

Put $r=4Dq$. Each degree in $T$ is at most $2Dq$. The determinant hence has size at most $r$, giving $\deg R\le Dr$ and $\|R\|_1\le p^{Br}$. For the norm bound use the triangle inequality on the determinant and bound the permanent of the matrix of entry norms by the product of its row sums. Each row sum is at most $p^B$. The degree-zero choice satisfies the same bounds.

For two support vectors $\nu,\mu$ of $R$ we now have
$$|(\nu-\mu)\cdot e|=|(\nu-\mu)\cdot(e-Kv/q)|\le 2DrK/(qQ)\le \epsilon K/2.$$
If $R(p^e)\ne0$, divide it by the lowest occurring power of $p$, obtaining an integer divisible by $W$ of absolute value at most $p^{Br+\epsilon K/2}$. This excludes $K\ge 8BDq_0/\epsilon$. Include smaller $K$, and $K=0$, by listing the bounded integer points.

If $R(p^e)=0$, order its monomials by $\nu\cdot e$, retaining ties. There are at least two. The largest exponent exceeds the second largest by at most $\log_p\|R\|_1\le Br$ since otherwise the leading integer term could not cancel. Thus
$$(\nu-\mu)\cdot e=c,\qquad 0\ne\nu-\mu,\qquad
\|\nu-\mu\|_\infty\le4D^2 q_0,\qquad |c|\le4BDq_0.$$
Listing all such integer equations finishes the construction. Its bounds can be large.

## Restricting to successive hyperplanes

Here is the iteration used below. Suppose a fixed collection of integer Laurent polynomials has no common nonunit factor upon restriction to any connected torus translate of codimension less than $k$ that contains a positive real point, including the whole torus. Suppose $1\le k\le n$. Exponents at which their common divisor away from $p$ exceeds $p^{\epsilon K}$ lie in an effective finite union of rational affine spaces of codimension at least $k$, simultaneously for all $p$.

On each affine space with integer points write $e=e_0+M w$, using an integer origin and a basis of its saturated lattice of differences. These are effective by integer linear algebra. The corresponding monomial map, translated by $p^{e_0}$, parametrizes a connected torus translate isomorphically because the lattice is saturated. An integer left inverse to $M$ gives $\|w\|_\infty\le A K$ for an effective $A$ whenever $K\ge1$. Clear Laurent monomials and powers of $p$ in the restricted polynomials. This preserves every common divisibility away from $p$. The degrees are bounded by a constant from $M$ and the original supports. The coefficient norms are bounded by $p^{B'}$ with $B'$ independent of $p$, since every coefficient is a sum of the original integer coefficients times powers of $p$ with bounded exponents.

For several polynomials, discard zero members. A Laurent unit gives a direct bound $p^{B'}$. Otherwise choose one nonconstant member $F$ and an integer linear combination $G$ of the others that avoids every nonmonomial irreducible factor of $F$. If the others are $G_0,\ldots,G_d$, try $\sum t^i G_i$ for integer $t$. Each factor excludes at most $d$ values since at least one $G_i$ is nonzero modulo that factor. Thus $t\le d\deg F$ suffices among the nonnegative integers, so even a base-dependent choice preserves uniform norm bounds. Apply the preceding hyperplane construction in $w$ with $\epsilon/A$. Its hyperplanes cut the current space properly, and bounded points may be kept as points. Repeating at most $k$ times gives the claim.

## Applying the bound to sparse rows

Use the [factor assignments and block geometry](sparse_prime_counts.md#the-factors-assigned-to-two-rows). Write $W=p^a$, $m=\sum_i p^{h_i}$, and the row avoiding $p$ as $(Wx,Wy,1+Wz)$ where $x+y+z=m$. Its assigned weight factor $U\ge3$ divides the row. Hence $x,y$ each use at least two units, and $z$ at least one. Split the assigned factors of $m$ into $S_{xy},S_x,S_y$ according to their two zero slots, with $S_{xy}$ on the two unshifted slots. Then
$$27m^2\le U^3(S_{xy}S_xS_y)^2,\qquad S_{xy}\mid D_0=\gcd(x,y,z),\qquad
U S_x S_y\mid 1+Wz.$$
Also $U\mid\gcd(x,y)$, $S_x\mid\gcd(x,y+z)$ and $S_y\mid\gcd(y,x+z)$. All these assigned factors are coprime to $p$.

For two disjoint sums with at least two units each and at least three in one block, the restriction hypothesis above holds for $k=2$. For three disjoint sums with at least two units each and at least three in one, it holds for $k=3$. The positive-coefficient statements here use exactly the [two-block and three-block character calculations](sparse_prime_counts.md#two-disjoint-sums). The common zero variety is an irreducible product of hyperplanes in a torus. Only a two-unit block fixes a nontrivial character, its coordinate ratio $-1$. The rank and sign exclude containment in a positive translate of the critical codimension. Normalizing one labelled coordinate to 1 preserves this calculation by homogeneity.

Take $\delta=1/10$. If $D_0>m^\delta$, the three-block bound places the inner pattern on codimension-three spaces, with two qualifications. A singleton gives $D_0=1$. If all three blocks are pairs, the [elementary three-gap argument](sparse_prime_counts.md#the-factors-assigned-to-two-rows) gives $D_0\le p^{H/3}+1$ in every feasible allocation. In detail, a zero gap gives $D_0\le2$. Otherwise $D_0>2$ divides $p^g+1$ for the gcd $g$ of the gaps, with all three gap quotients odd. Equal gaps would give $U\mid 1+p^g\mid m$, impossible. Unequal gaps force $3g\le H$.

In every case outside the codimension-three spaces, at least one of $U,S_x,S_y$ exceeds $m^\delta$. Otherwise the size bound would give $27m^2\le m^{7\delta}\max(m^{2\delta},4m^{2/3})$, impossible. A large $S_x$ or $S_y$ gives two disjoint sums of lengths at least two and three. A large $U$ gives the same when one of $x,y$ has at least three terms. The effective $k=2$ bound places all such inner patterns in codimension-two spaces since $m\ge p^H$.

If $U$ is large and $x,y$ are pairs, let $d>0$ be the gap of $x$. The congruence $p^d=-1\bmod U$ permits replacing $a$ by its least nonnegative residue $a'$ modulo $2d$. The three sums $p^{a'}x,p^{a'}y,1+p^{a'}z$ are still multiples of $U$. Their lengths are $2,2,s-4$, and the exponents above the constant coordinate are at most $3H$. Apply $k=3$ with $\epsilon=\delta/4$ in the $s-1$ coordinates $a'+h_i$. Subtracting the first coordinate and deleting it maps each resulting space to dimension at most $s-4$ in the inner space of dimension $s-2$. This proves the inner position assertion. If an integer dilation $b h$ belongs to one of these spaces, either the space is homogeneous or one defining integral equation has nonzero right side $c$. In the latter case $b\mid c$, so only finitely many dilations need be applied inversely to the space. This proves the dilation assertion with the same codimension.

For the full count per prime, $K_0=m\gcd(x,y)\le m^2$ is divisible by whichever of $U,S_x,S_y$ is large. Hence $\gcd(K_0,1+p^az)>p^{\delta H}$. For $H\ge2\log_2(s-1)$ we have $K_0\le p^{3H}$. The [elementary lcm spacing estimate](sparse_prime_counts.md#how-many-outer-shifts-can-supply-a-large-divisor) gives $O(1+L/H)$ possible shifts for an allocation at a fixed inner pattern and prime, with absolute effective constants for these choices of 3 and $\delta$. In a dyadic range $R\le H<2R$ the codimension-two spaces contribute $O_s(R^{s-4}(1+L/R))$ tuples. Summing gives $O_s((1+L)^{s-4})$. The codimension-three inner spaces permit $O_s((1+L)^{s-5})$ patterns and at most $L$ shifts. Bounded $H$ contributes $O_s(1+L)$. All constants are effective and independent of $p$.

## Combining primes when the shift is at most a fixed multiple of the span

Use the original $s-1$ exponent coordinates, writing $X,Y,Z$ for sums over the three disjoint labelled blocks, so their numerical values are $Wx,Wy,Wz$. For $U,S_x,S_y$ respectively, the block sums modulo the assigned factor equal the following triples of constants:
$$(0,0,-1),\qquad (0,1,-1),\qquad (1,0,-1).$$
For $S_{xy}$ they equal $(0,0,0)$. In each instance subtract these constants to get three polynomials. Each block equation $\sum t_i=c$ cuts out an irreducible hyperplane intersected with the torus, or has empty zero set. A block with $c=0$ fixes a character only if it has two terms, when the ratio is $-1$. A block with $c\ne0$ fixes a character only if it is a singleton, when the coordinate equals $c$. To see the latter assertion for two or more terms, eliminate the last variable. The other coordinate variables and $c$ minus their sum are distinct nonassociate irreducibles.

For $U$, at most two blocks fix characters. Three would require $x,y$ both pairs and $z$ a singleton, using only five positive units. For $S_x,S_y$ the block with constant $1$ has at least two terms and fixes no character. In these cases every fixed generator has value $-1$. For $S_{xy}$ the same statements hold except for three pairs, already bounded separately by $S_{xy}\le 2m^{1/3}$. A singleton with constant zero gives a unit polynomial and bounds the assigned factor by 1.

In the nonempty cases the product zero variety has codimension three and is irreducible. On a translate of codimension zero or one, a common nonunit factor would give a component of dimension greater than that of this variety. On a translate of codimension two it would force the whole variety into the translate by equality of dimensions. A connected codimension-two translate fixes a saturated rank-two lattice of characters. If the variety fixes two generators, they are supported on disjoint blocks and generate a saturated lattice. Containment would force equality of the two lattices, hence a negative character value on a translate containing a positive point. Thus the restriction hypothesis holds for $k=3$.

When $a\le CH$, the largest exponent is at most $(1+C)H$. Apply the effective $k=3$ result with $\epsilon=1/(10(1+C))$. Outside the resulting union of spaces, all assigned factors are at most $m^{1/10}$, except for the separately bounded three-pair factor. The same size contradiction proves the combined assertion.

## An explicit family near equally spaced positions

For the positions $h=(0,b(d^n+1),b(2d^n+d),\ldots,b(nd^n+d^{n-1}))$ in the introduction, put $P=p^b$. Write $m=\sum p^{h_i}$. Three numerical gcds away from $p$ bound the assigned divisors. Namely $U$ divides the first of
$$\gcd(x,y),\qquad \gcd(x,m-x),\qquad \gcd(y,m-y),$$
and $S=S_{xy}S_x S_y$ divides the product of the last two. Indeed at every power in $S$, at least one of $x,y$ vanishes modulo that power and the power divides $m$. If all three gcds away from $p$ are at most $L_0$, (1) of the [two-row size argument](exponent_fan.md#the-two-divisors-to-bound), with the assigned $p$-power of the predecessor canceled, requires $27m^2\le L_0^7$.

For each of these gcds use the two disjoint block sums as the linear polynomials $F,G$, with the lowest labelled position represented by 1. They are coprime as Laurent polynomials by disjoint supports in independent variables. Use the fixed projection direction $v=(1,2,\ldots,n)$. Its resultant $R$ has degree at most $2n$ and norm at most $(n+1)^{2n}$ because each degree in $T$ is at most $n$. All monomials have the same $v$-weight, so the common factor contributed by $d^n v$ in the exponents can be removed. Only the residual exponents $(1,d,\ldots,d^{n-1})$ remain in the base $P$.

Put $d_0=2n(1+\lceil\log_2(n+1)\rceil)+1$ and take $d\ge d_0$. For two different monomials of $R$, the highest differing position $i$ is at least two since their $v$-weights agree. Orient their difference to be positive at $i$. The negative coefficients below $i$ have total absolute value at most $2n$. Therefore their residual exponents differ by at least $(d-2n)d^{i-2}\ge d-2n>\log_2\|R\|_1$. The largest term cannot cancel at $P\ge2$. The whole span of its exponents is at most $2n d^{n-1}$, so each common divisor is bounded by
$$L_0=(n+1)^{2n}P^{2n d^{n-1}}.$$
Yet $m\ge P^{n d^n}$ and $d^{n-1}(d-7)>7\log_2(n+1)$ since $n\ge5$ and $d\ge d_0$. Thus $27m^2>L_0^7$.

**Seven units at every integer $d\ge2$.** For $s=7$ one has $n=5$ and $d_0=41$. For each $2\le d\le40$, the second verifier below treats the 105 possible pairs of unshifted blocks. They use two and two units, or two and three, with at least one unit left for $z$. For each of the three gcds it collects equal powers in the projected univariate polynomial, strips the lowest power, and records the degree $k_i$ and coefficient norm $c_i$. The polynomial has no integer root $P\ge2$. This follows either by domination of the leading term or by testing the positive divisors of its nonzero constant coefficient. With $J=5d^5+d^4$ and $(w_1,w_2,w_3)=(3,2,2)$ the verifier checks
$$g=2J-\sum_i w_i k_i\ge0,\qquad 27\cdot 2^g>\prod_i c_i^{w_i}.$$
The same strict inequality holds with $P$ in place of 2. As $m\ge P^J$, every allocation violates the necessary size bound. This checks 4,095 allocation and parameter pairs using 5,460 distinct polynomial specializations, with $a,b$ and the prime unrestricted throughout.

## An exact projection to reproduce

For $F=1+X$, $G=1+Y+Z$ and $v=(3,2,4)$, the resultant, up to sign, is
$$X^4+X^2Y^3+Z^3-3X^2YZ.$$
All terms have weight 12. At $(X,Y,Z)=(p^a,p^b,p^c)$ the expression vanishes exactly at $4a=2a+3b=3c$, by arithmetic-geometric mean. On this ray, with $(a,b,c)=(3j,2j,4j)$ and $t=p^j$, the numerical gcd is exactly $t^2-t+1$, by factorization. Off the ray the common divisor divides the nonzero integer obtained by stripping the smallest power of $p$ from the resultant.

Run the following from the repository root with SymPy installed. The first program checks 12,960 signed-exponent specializations including 251 vanishing projections, as well as 1,000 exact gcds along the ray. The second verifies the explicit seven-unit family. The third checks the period and shift arithmetic used in the counting theorem.

~~~sh
python3 check_sparse_projection.py
python3 check_perturbed_positions.py
python3 check_sparse_prime_counts.py
~~~
