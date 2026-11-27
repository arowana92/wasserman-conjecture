# Binary orders of possible counterexamples with a bounded denominator

Consider a possible counterexample with a nonendpoint row ordered as $(x,y,z)$. Put $W=\gcd(x,y,z-1)$ and $d=(N-1)/W$.

**Theorem.** For $d\le1000$, at most $O(H^{0.778})$ of any $H$ consecutive binary orders contain such a counterexample. Orders can be taken on $N$ or on $W$. If $W$ is a cube the bound is $O(H^{0.721})$. Both bounds are uniform in the first order and allow arbitrary prime factors of $W$.

In particular the first bound counts orders containing any failing weight of the form $1+dq^e$, simultaneously over $1\le d\le1000$, all primes $q$, and all $e\ge1$. A window of $H$ orders contains a run of $\gg H^{0.222}$ consecutive orders free of these weights, for sufficiently large $H$. Taking the binary orders of the **base prime** instead, at most $O(H^{0.9558})$ orders in each window of length $H$ contain any failure at $d\le1000$ and $1\le e\le H^{1/10}$.

For the separate family of even weights $7t+1$ with $t=p^a q^b$, $p,q$ prime and $a,b\ge1$, the number of possible exceptional binary orders of $t$ in any such window is
$$
 O\big(H^{\log_2\lambda}(1+\log H)^5\big),\qquad
 \lambda^3=\lambda^2+1,\quad \lambda>1 .
$$

## Two independent bases

Suppose $W$ is accepted by digit graphs in two multiplicatively independent integer bases $p,r$. Suppose accepted expansions allow padding by high zeros, and the numbers of possible length-$k$ internal words are at most $C_p p^{sk}$ and $C_r r^{tk}$ at the respective bases. Read each word in the order of significance when forming its integer value. If $s+t<1$, then for every $\epsilon>0$, the number of occupied base-$p$ orders in any $H$ consecutive orders is $O_\epsilon(H^{s+t+\epsilon})$.

Choose an integer block length $L$. For $b=p,r$, let the alphabet in base $B_b=b^L$ consist of the length-$L$ admissible words together with zero. Let $K_b$ be the compact subset of $[0,1]$ whose digits belong to this alphabet, allowing either expansion at endpoints. Its dimension is $s_b=\log_{B_b}|\mathcal D_b|$. For any given $\eta>0$, $L$ can be chosen so that $S=s_p+s_r\le s+t+\eta<1$. Put $n_b=1+\lfloor\log_{B_b} W\rfloor$ and $u=W/B_p^{n_p}$, $v=W/B_r^{n_r}$. Padding the top of the integer to a block boundary and appending zero fractional blocks shows that
$$
 (u,v)\in K_p\times K_r,\quad u\ge 1/B_p,\quad v\ge 1/B_r,\qquad
 n_p\alpha=\log_{B_r}(v/u)\pmod1,\quad \alpha=\log_r p .
$$

We use a uniform strip estimate for this product. For a line of slope in a fixed compact subinterval of $(0,\infty)$, the part of the product within distance $h$ of the line meets $O_\xi(h^{-\xi})$ square grid cells of side $h$, for every $\xi>0$ and $0<h\le1$. Indeed, the natural product probability gives every ball of radius $h$ centered on its support mass $\gg h^S$. Its projection in the normal direction gives every interval of length $O(h)$ mass $O_\xi(h^{S-\xi})$, uniformly for the specified slopes. These are the regularity and the uniform projected interval bounds proved for restricted-digit products by Corso and Shmerkin. The latter uses multiplicative independence of $B_p,B_r$. Pack disjoint balls of radius $h/3$ with centers in the strip, enlarge the strip by $h$, and divide the upper mass bound by the lower one. A maximal packing covers the required points with balls of radius $h$, proving the cell estimate. [Corso and Shmerkin, proof of Theorem 1.15, equation (6.3) and the projected dyadic interval estimate following it](https://arxiv.org/pdf/2409.04608)

Count $H$ consecutive indices $n_p$, starting at $n_0$. Take a convergent $a/q$ to $\alpha$ with $q\le H<q'$ for its next denominator. Write $\sigma=|q\alpha-a|<1/q'$ and $R=1/q$. Distinct integer bases that are multiplicatively independent satisfy $\sigma\gg q^{-K}$ for a constant $K$ depending only on the bases, by a lower bound for a nonzero linear form in two logarithms. The few small denominators can be absorbed in the constants. [Gouillon, Corollary 2.3](https://www.numdam.org/article/JTNB_2006__18_1_125_0.pdf)

Group the indices by $n_p-n_0$ modulo $q$. All phases of a group are within $R$ of its corresponding point of the translated $q$-grid, since $H|\alpha-a/q|<R$. They advance by $\pm\sigma$ and lie on an arc of length at most $(H/q)\sigma<R$. Set $w=\min(R,(1+H/q)\sigma)$, so $\sigma\le w\le R$. Cover the product by its digit rectangles with sides between $R/B_b$ and $R$, using the first cylinder depth reaching scale $R$ in each coordinate. There are $O(R^{-S})$ rectangles. The logarithmic ratio is Lipschitz on the positive rectangle of admissible mantissas, so each digit rectangle contributes to only $O(1)$ groups by grid spacing.

For a given group, pull back its arc of length at most $w$ under the logarithmic ratio. Within the admissible mantissa rectangle this is contained in a bounded number of strips of width $O(w)$ about lines $v=\beta u$, with $\beta$ bounded above and away from zero. The bounded number accounts for reduction modulo one. In a fixed scale-$R$ digit rectangle, rescale each coordinate by its cylinder length. The product there is a translated diagonal copy of $K_p\times K_r$ before rescaling, the new slope remains in a fixed compact subinterval of $(0,\infty)$, and the new strip width is $O(w/R)$. The strip estimate leaves $O_\xi((R/w)^\xi)$ product cylinders of scale comparable to $w$ inside this rectangle. A grid cell of comparable size meets only a bounded number of these cylinders per coordinate.

Subdivide each surviving cylinder pair to scale $\sigma$. This costs $O((w/\sigma)^S)$, since each extra full block has precisely $|\mathcal D_b|$ choices. Each fine rectangle projects to an arc of length $O(\sigma)$ and thus accounts for $O(1)$ indices in this group by the progression spacing (split at the circle boundary if needed). The resulting count is at most
$$
 O_\xi\big(q^S(R/w)^\xi(1+H/q)^S\big)
 \ll_\xi H^{S+K\xi}.
$$
Choosing $\eta,\xi$ sufficiently small gives the claimed exponent. An interval of base-$p$ orders meets at most $H+1$ block orders and each occupied block covers at most $L$ base-$p$ orders, a constant for fixed $\epsilon$.

## Carry graphs and the sixteen rows

The [primitive denominator theorem](shifted_directions.md) and the [exclusions at 364 and 880](denominator_residue_hulls.md) give the rows below. Each indicated pair of primes requires the first row $(aW,bW,cW+1)$ to add without carries at both bases. The same reduction forces a nonendpoint row with exactly $W=q^e$ in any failing weight $1+dq^e$ in the stated cofactor range. Only the first two rows at $986$ are possible for cube scales. [Cube congruences](common_power_exponents.md)

| $d$ | $(a,b,c)$ | Bases | Upper exponents in units of $1/10000$ |
|---:|---|---|---|
| 380 | 130,247,3 | 5,13 | 2848,3875 |
| 589 | 133,434,22 | 7,19 | 2564,4238 |
| 629 | 259,357,13 | 7,17 | 2529,3713 |
| 646 | 133,476,37 | 7,17 | 2983,3965 |
| 703 | 133,518,52 | 7,19 | 2457,4255 |
| 710 | 175,497,38 | 5,7 | 1591,2393 |
| 738 | 117,533,88 | 3,13 | 1987,3702 |
| 780 | 224,455,101 | 2,5 | 1880,2068 |
| 850 | 175,476,199 | 5,7 | 1622,2593 |
| 874 | 209,253,412 | 11,19 | 3360,4416 |
| 946 | 301,539,106 | 7,11 | 2156,3429 |
| 950 | 275,627,48 | 5,11 | 2405,3520 |
| 962 | 143,814,5 | 13,11 | 3328,3495 |
| 986 | 119,203,664 | 7,17 | 2839,4364 |
| 986 | 203,714,69 | 7,17 | 2473,4093 |
| 986 | 377,442,167 | 13,17 | 3481,4235 |

Delete states which cannot reach acceptance from the carry graph, and let $M$ count the labeled transitions between remaining states. For $v=M^{300}\mathbf1$, the checker verifies $10^6Mv\le Fv$ and $F^{10000}\le(10^6)^{10000}b^j$ for the printed exponent numerator $j$ in base $b$. Iteration bounds all internal words by $C b^{jk/10000}$, summing over starting states and absorbing lengths below 300. An accepting state permits more high zeros because its remaining carries already add freely. The largest exponent sum is $7776/10000$, and at cubes it is $7203/10000$. Applying the theorem and summing the sixteen bounds gives the result. For each fixed $d$ a base-$p$ order of $W$ intersects only $O(1)$ binary orders of $W$ or of $1+dW$, and a window of $H$ binary orders meets $O(H)$ such base-$p$ orders. A window containing $M$ exceptional orders leaves a run of at least $(H-M)/(M+1)$ good orders. For base-prime orders and a fixed $e$, the window of orders of $q^e$ has length $eH$; division of an occupied order by $e$ and rounding down recovers the order of $q$. Summing $(eH)^{0.778}$ for $e\le H^{1/10}$ gives the exponent $0.9558$.

## Even weights at cofactor seven

In this family, $6t$ and $t+1$ must add without carries in binary. The base-seven digits of $t$ are all nonzero and weakly increasing from the most significant end. There are at most $\binom{k+5}{5}$ possible leading words of length $k$. [Factor-seven reduction and digit inequality](predecessor_five_seven.md)

The binary transition matrix for the multipliers $(6,1)$ and offsets $(0,1)$ has characteristic polynomial
$$ x^3(x-1)(x^3-x-1)(x^3-x^2-1). $$
The two cubics are characteristic polynomials of primitive nonnegative companion matrices. Their Perron roots are respectively smaller than $\lambda$ and equal to $\lambda$ (evaluate the first cubic at $\lambda$). The unique eigenvalue of largest modulus is simple. Thus the number of any length-$k$ internal binary words is $O(\lambda^k)$. Write $\delta=\log_2\lambda$. In normalized binary mantissas bounded away from zero, every interval of length $R$ meets $O(1)$ binary cylinders of comparable width, each with $O((R/h)^\delta)$ admissible descendants at comparable width $h\le R$. This still holds in the logarithm of the mantissa because its derivative is bounded above and away from zero.

Apply the grid argument with $\alpha=\log_7 2$ and choose $k$ so that $7^{-k}\asymp\sigma$. The bound $\sigma\gg H^{-K}$ gives $k=O(1+\log H)$. Orders with too few base-seven digits contribute only $O(1+\log H)$ initial binary orders. Fix one of the $O((1+\log H)^5)$ eligible base-seven prefixes. Up to an arc error $O(\sigma)$ the phase then belongs to a fixed translate of the negative logarithm of the binary mantissa set. Finite binary inputs can be extended at the fractional end with zeros; an internal word crossing that end has a prefix in the original language and a suffix of zeros. Summing the geometric bounds over prefix lengths preserves $O(\lambda^k)$. The global cover at scale $R$ leaves $O(q^\delta)$ groups. The local cover within each group at scales $w,\sigma$ leaves $O((1+H/q)^\delta)$ indices, including the phase error by enlarging arcs. Sum over the seven-adic prefixes to obtain the stated bound.

Run the certificate with Python, python-flint and SymPy.

~~~~sh
python3 check_two_scale_orders.py
~~~~
