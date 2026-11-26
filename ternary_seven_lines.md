# One-parameter families for the seven-unit case at odd primes

Let an odd prime $p$ divide $N-1$ and suppose $N$ has base-$p$ digit sum seven. Expand the six positive-position units with multiplicity as
$$N=1+\sum_{j=1}^6 p^{a+h_j},\qquad a\ge1,\quad 0=h_1\le\cdots\le h_6=H.$$
Across all these primes, the normalized patterns $(h_1,\ldots,h_6)$ that could give a counterexample at any shift lie on **finitely many rational affine lines**, together with a finite set. In particular at most $O(1+B)$ patterns of span $H\le B$ can fail. Counting the actual weights through $X$, including the outer shift, gives
$$O\big((1+\log X)(1+\log(1+\log X))\big).$$
The implied constants are absolute. This strengthens the [quadratic logarithmic bound](odd_seven_count.md). The large-span assertion below is at $p=3$. Only finitely many patterns for all $p\ge5$ combined remain by the [large-span theorems](seven_small_prime_span.md).

## Two factors must both be large

Assume $p=3$, $H\ge10^{24}$ and an exception. The [seven-unit reduction](seven_small_prime_span.md#the-three-common-factors) gives three inner blocks
$$ x=X^i(1+X^d),\quad y=X^j(1+X^e+X^f),\quad z=X^k,\qquad d>0,\quad 0\le e\le f. $$
Here and below an integer divisor of a polynomial value refers to evaluation at 3. For assigned positive odd integers prime to 3 we have
$$ U\mid (1+3^d),y(3),\quad A\mid (1+3^d),(y+z)(3),\quad C\mid y(3),(x+z)(3), $$
$$ U\ge7,\qquad \gcd(C,3^{3d}-1)=1,\qquad 27\cdot3^{2H}\le U^3 A^2 C^2. \tag{1} $$
We prove $A,C\ge3^{H/1000}$. All logarithms below are to base 3.

### Simultaneous bounds on the divisors of the binomial

List the four exponents of $y+z$ in increasing order after subtracting the least, as $a_1=0,a_2,a_3,a_4$. Let $H_0=\max(d,a_4)\le H$. For $H_0\ge16$ put $Q=\lfloor H_0^{1/4}\rfloor$. Apply the [four-coordinate approximation](binomial_sparse_gcd.md#the-remaining-integer-factor) to $(d,a_2,a_3,a_4)$, obtaining a denominator $q\le Q^3$, primitive nonnegative slopes $(n,v_2,v_3,v_4)$ with $v_1=0$ and a scale $\lambda$. Their maximum $h\le q$ satisfies $\lambda h=H_0$. Every error $|d-\lambda n|$, $|a_j-\lambda v_j|$ is at most $H_0/(qQ)$. Round nearest with a common tie rule, so the four slopes are ordered.

Use $f_0=1+3^d T^n$ and the two polynomials formed from $3^{a_j}T^{v_j}$, one on the three labels of $y$, the other on all four. Strip lowest powers of $T$. Write $d_U,d_A$ for their respective gcd degrees with $f_0$, and $D=3d_U+2d_A$. Then
$$ 3\log U+2\log A \le (D/h)H_0+50H_0^{3/4}. \tag{2} $$
Indeed the coefficient-minor bound gives scalar loss at most $8^q3^{3H_0/Q}$ for each gcd. The three-label sublist has at most the residual spread of the four-label list, even if it lacks the original constant. Stripping a monomial adds a fixed weight to each column and preserves its residual spread. For $n>0$ a primitive gcd of degree $m$ divides the primitive binomial with constant 1. Each root has modulus $3^{-d/n}$, so its value at 1 is at most $2^m 3^{md/n}$. Since $D\le5n$, the error replacing $D d/n$ by $D\lambda$ is at most $5H_0/(qQ)$. The combined error is bounded by $20q+20H_0/Q\le50H_0^{3/4}$. For $n=0$ the scalar bound suffices (also $d\le H_0/(qQ)$). For $H_0<16$, $U,A\le1+3^{15}$ already force the desired lower bound on $C$.

At a root $T=\xi3^{-d/n}$ put $b_j=a_j-dv_j/n$, so the four terms have phases $\xi^{v_j}$ and magnitudes $3^{b_j}$, with $\xi^n=-1$. We may exclude $d_A=n$. In that case Gauss's lemma on the original integer polynomials gives $1+3^d\mid (y+z)(3)$ by evaluation at $T=1$, contradicting $U\mid y(3)$ and $U\nmid z(3)$.

Count first the quad roots which cancel within each group of equal height. They give a matching of opposite phases with equal heights on each edge. For a fixed matching with slope gaps $w_1,w_2$, its roots, when nonempty, solve $\xi^{\gcd(n,w_1,w_2)}=-1$, as in the [opposite-matching proof](seven_large_bases.md#common-roots-at-the-unshifted-pair). A zero gap cannot be opposite. At most two matchings are admissible by the residue classes modulo $2^{v_2(n)+1}$. A full-binomial matching was excluded above. Hence the number $d_{A,0}$ of these roots is at most $7n/15$.

If $d_U>0$, put $g=\gcd(n,\{v_j-v_l:j,l\text{ in }y\})$. After subtracting a common slope, the [three-term phase classification](binomial_trinomial_gcd.md#locating-every-common-root) gives $d_U\le2g$. For equal heights on $y$, the normalized slopes, starting from the smallest, have residues $0,2,4\bmod 6$ in some order, and $n/g$ is odd. Thus $4g\le h$, and no pair inside $y$ has the 2-adic valuation of its slope difference equal to $v_2(n)$. A matching would require such a pair. For unequal heights on $y$, the two lowest phases have ratio a primitive sixth root, and the highest has phase of order 12 relative to either lowest. Thus $6\mid n/g$ since the normalized root has order divisible by 12 and its $(n/g)$th power is $-1$. An equal-height matching must pair the two lowest inside $y$. For their slope gap $w$, the order of $\xi^w$ is divisible by 3, so $v_3(w)<v_3(n)$ for the root of the **trio**. The number of **quad** roots making that same pair opposite is at most $\gcd(n,w)\le n/3$. Consequently $3d_U+2d_{A,0}\le5h/3$ in all cases, including $d_U=0$.

All remaining quad roots give minimal mixed-height relations. For any fixed heights the [ternary four-term table](quadrinomial_small_bases.md) gives either order 2 or order 12 for the group of relative phases, with just one order for those heights. As $v_1=0$, adjoining $\xi^n=-1$ preserves that group; primitivity shows it contains $\xi$. There are at most four such roots. Thus $D\le5h/3+8$. For $25\le h\le48$ this gives $D<2h$ and integrality saves at least $1/h$ in $D/h$. For $h>48$ the ratio is at most $11/6$. Unless a mixed root occurs with $h\le24$, we have $D/h\le2-1/48$, which by (1) and (2) gives $C>3^{H/1000}$.

### The two cases of equality

The [exact verifier](check_ternary_seven_lines.py) enumerates the remaining directions. For $1\le n\le24$ take all $0=v_1\le v_2\le v_3\le v_4\le24$ primitive with $n$. In the mixed case the tuple $\delta_j=2(b_j-b_1)$ belongs to the permutations of the four height patterns in the ternary table after reanchoring. Keeping only phases from that table with $\xi^n=-1$ covers every mixed root. The integers $2d,(2a_j-\delta_j)_j$ are proportional to the primitive slopes, so $B=d/n$ is half-integral. Compute the two gcd degrees with $B$ replaced by $b/2$, $b=0,1$ of the same parity as $2B$. This just rescales $T$ by a power of 3. Retain integral powers and discard any three coincident terms (equal slopes and equal heights would occupy one ternary digit). After excluding 1,432 full-binomial configurations, there are 64,008 labelled allocations. The only ones with $D\ge2h$ are

| $n$ | four slopes | doubled heights | singleton slope | $d_U,d_A$ |
|---:|---|---|---:|---|
| 3 | $0,1,2,4$ | $0,2,0,0$ | 1 | $2,1$ |
| 3 | $0,2,3,4$ | $0,0,2,0$ | 3 | $2,1$ |

All other small directions have $D/h\le2-1/24$. In the table $B$ is a positive integer. Put $Z=3^B$ and $\Phi_3(Z)=Z^2+Z+1$, $\Phi_6(Z)=Z^2-Z+1$. Up to a power of 3 the trio equals $1+Z^2+Z^4=\Phi_3(Z)\Phi_6(Z)$, and the binomial is $1+Z^3=(1+Z)\Phi_6(Z)$. Since $\Phi_3(-1)=1$, $U\mid\Phi_6(Z)$. Any prime of $A$ dividing $\Phi_6(Z)$ would divide $y(3)$ and $(y+z)(3)$, so $A\mid1+Z$. Likewise $C$ is coprime to $1+3^d$ by $C\mid (x+z)(3)$; it is coprime to $\Phi_3(Z)$ by (1) and $\Phi_3(Z)\mid Z^9-1$. Hence $C=1$. Now $U^3A^2<4Z^8$ while $H\ge4B$, contradicting (1). This finishes the lower bound on $C$.

### The divisor belonging to the four-term block

Suppose $A<3^{H/1000}$. The [two-thirds bound](two_trinomial_cuts.md) gives $C\le4\cdot3^{2H/3}$. Therefore $\log U>0.22H$ by (1). Apply the [single binomial–trinomial approximation](binomial_trinomial_gcd.md#bounding-the-remaining-scalar) to $d,e,f$. The bounded internal heights are immediate by $U\le1+3^d$. The line gcd has positive degree $\mu\le2$ and there are primitive nonnegative slopes $n,u,v$ with $n>0$, $u\le v$, and a positive half-integer $B$ such that
$$ d=Bn,\quad e=Bu+s/2,\quad f=Bv+t/2,\qquad
(0,s,t)\in\{(0,0,0),(0,1,0),(0,0,1),(0,-1,-1)\}. $$
Here $\log U\le\mu B+6H^{2/3}+2$ by the root modulus argument used in (2). In particular $B>H/9.2$, $h_0=\max(n,v)\le9$ and $\mu/h_0\ge0.21$ since $Bh_0\le H+1$. Enumerating all these primitive slopes, offsets and scale parities as above gives precisely the following list. All degrees are two.

| $n,u,v$ | $s,t$ | $2B\bmod2$ |
|---|---|---:|
| $3,2,4$ or $3,4,8$ | $0,0$ | 0 |
| $6,1,2$ | $1,0$ | 1 |
| $6,2,7$ | $0,1$ | 1 |
| $6,5,7$ | $-1,-1$ | 1 |
| $9,2,4$ or $9,4,8$ | $0,0$ | 0 |

The first profile $(3,2,4)$ forces $C=1$ by exactly the two coprimalities for $\Phi_3,\Phi_6$ above, independent of the singleton. All others have $h_0\ge6$. Also $\log C>0.24H$ follows from (1), the assumed small $A$ and the [half-span bound](odd_seven_count.md#a-half-span-bound-on-the-first-common-divisor) on $U$.

Put $F(S)=1+S^n$, $Y(S)=1+3^{s/2}S^u+3^{t/2}S^v$, $\tau=i-k$. Choose coprime $1\le q\le32,r\in\mathbb Z$ with $|q\tau/B-r|\le1/32$; then $|r|\le10q$. Write $\eta=\tau-Br/q$. At $T=3^{B/q}$ the divisor $C$ divides the two values of
$$Y(T^q),\qquad 1+3^\eta T^r F(T^q) \tag{3}$$
in the ring with 3 inverted. A nonzero resultant after clearing monomials and denominator powers would give $\log C\le 600q+(H+2)/32$ by the [two-variable resultant estimate](odd_seven_count.md#a-second-restriction-on-the-six-positions) (the two degrees are at most $9q,19q$). This is too small. The common-root height estimate there gives $|\eta|\le7(n+|r|/q)+1\le140$.

Write $d_C$ for the common degree in (3). It is at most $v-\mu$. Indeed the $\mu$ distinct roots shared by $Y,F$ cannot satisfy the second equation. For each remaining root $S$ of $Y$, (3) fixes both $T^q=S$ and $T^r$, admitting at most one $T$ by coprimality. All roots considered are nonzero, and multiplicity from $Y(T^q)$ equals that from $Y(S)$. The [coefficient-minor and norm bound](sharp_trinomial_gcd.md#bounding-the-numerical-gcd) gives
$$ \log C\le d_C B/q+50000. \tag{4} $$
In detail, the two inputs have coefficient exponent spreads at most 1 and 140, before collecting like terms. Clear Laurent monomials and the least coefficient powers. At each embedding, each coefficient minor is bounded by $3^{(9+19)q+19q+9\cdot140q}$ using the column sums. There are at most $v+1\le10$ coefficients after eliminating to degree $d_C$. The resulting gcd multiple evaluates to a nonzero algebraic integer, divisible by $C$ away from 3. The norm gives (4).

For $q\ge2$ we have $3\mu+2d_C/q\le2\mu+v\le5h_0/3$. It follows that $3\log U+2\log C\le(5/3)(H+1)+18H^{2/3}+100010$, contradicting (1) with the assumed small $A$.

For $q=1$ let $b_0\in\{0,1/2\}$ be the fractional part of $B$. Write $T=3^{b_0}V$ and $J=\eta+b_0r=\tau-(B-b_0)r$. This is an integer of absolute value at most 160. The two transformed polynomials have rational coefficients. For all seven profiles, $-11\le r\le11$, $-160\le J\le160$, the verifier computes the degree $d_C$ of their gcd after clearing denominator powers and Laurent monomials. All 51,681 checks give
$$3\mu+2d_C\le (3/2)h_1,\qquad h_1=\max\big(h_0,\max(0,r+n)-\min(0,r)\big).$$
Here $Bh_1\le H+281$: the actual positions of $x+z$ relative to $k$ are $0,\tau,\tau+d$, within $|\eta|$ of $B$ times $0,r,r+n$. Combining this estimate with (4) again contradicts (1). This proves the claimed lower bound on $A$.

## One scale for all six positive units

The [common-scale theorem](sharp_trinomial_gcd.md#consequence-for-seven-unit-counterexamples) gives $L\ge H/600000$ with
$$ 2d=Ln,\qquad 2e=Lu+s,\qquad 2f=Lv+t,\qquad |s|,|t|\le2. $$
Use this scale here regardless of the size of $U$. Put $K=600000$ and $D_0=2K+1$. All three slopes are at most $D_0$. Apply the [translation proof](odd_seven_count.md#a-second-restriction-on-the-six-positions) twice, to $k-j$ using $A$ and to $i-k$ using $C$. Use Dirichlet denominators up to 2000 in place of 32. Each nonzero resultant would bound the associated logarithm by $50(K+1)q+(H+2)/2000<H/1000$. Thus both pairs have a common root. For **each** translation $w$ there are integers $1\le q\le2000$, $|r|\le D_0q$ with
$$|2qw-rL|\le2\big(7(nq+|r|)+q\big)<10^{11}.$$
The offsets on the left are bounded integers, independently for the two translations. Fix these integers and the finitely many choices of all the displayed slopes, small offsets and denominators. Then $d,e,f,k-j,i-k$ are rational affine functions of the single integer $L$. Choosing which block starts at the minimum position fixes the common translation; a finite number of orderings gives the sorted pattern. This proves the affine-line assertion. There are $O(1+B)$ possible normalized patterns for span at most $B$, since $L\le2d\le2B$.

For $B=\log_2 X$ and dyadic spans $R\le H<2R$, there are $O(R)$ patterns. The [shift-spacing theorem](sparse_prime_counts.md#an-explicit-shift-bound-for-seven-units) gives at most 315 possible exceptional shifts in any interval of length $H/9$. Thus each bin contributes $O(R+B)$ weights; summing and including bounded spans gives the stated bound. The conjecture for the remaining affine families, as well as the general three-row case, stays open.

## A direct gap test

For any odd predecessor prime with seven digit units and $H\ge10^{24}$, **any** pair of positive-exponent positions at distance $g$ satisfying
$$20{,}000{,}000<g<H/10^{11}$$
proves the conjecture at that weight for every outer shift. Counterexamples at these spans would require $p=3$. In the preceding relations one of the two denominators can be chosen at most 32, since the [half-span argument](odd_seven_count.md#a-second-restriction-on-the-six-positions) gives $A\ge3^{H/9}$ or $C\ge3^{H/9}$. Use that choice for the larger divisor and the bound 2000 for the other one.

The difference of any two positions is thus of the form $cL+b$ with $2q_Aq_Cc$ integral. Its constant term has $|b|<20{,}000{,}000$: a difference of block starts uses each translation at most once, with absolute residual at most 6 for $A$ and $7(n+|r|/q)+1\le14D_0+1$ for $C$. The internal offsets add at most two. If $c=0$ the gap is small. Otherwise its size is at least $H/(600000\cdot2\cdot32\cdot2000)-20{,}000{,}000>H/10^{11}$. This proves the test.

For example at every odd prime $p$, every $a\ge1$ and every integer $t\ge5\cdot10^{11}$, it proves the conjecture for
$$ N=1+p^a\big(1+p^t+p^{2t^2}+p^{3t^2}+p^{4t^2}+p^{5t^2}\big). $$
Here $H=5t^2$ and $g=t$. More generally, any sequence of such patterns with one pair of positions at a distance $g$ satisfying $g\to\infty$ and $g/H\to0$ is covered eventually.

The first two verifiers use python-flint. A separate standard-library verifier counts binomial roots in split finite fields and repeats the translation bounds using the polynomial Euclidean algorithm modulo a prime. Run all three from the repository root.

~~~sh
python3 check_quadrinomial_small_bases.py
python3 check_ternary_seven_lines.py
python3 check_ternary_seven_lines_cross.py
~~~
