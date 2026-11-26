# A two-thirds bound for two trinomials and a reduction at bases 2 and 3

Let $p$ be an odd prime and let $P,Q$ each be a sum of three monomials with coefficient one, counting repetitions. Independently remove the lowest power of $X$ from each, and suppose the resulting polynomials are different. Let $H$ be the greater of their two degrees. Then
$$\gcd(P(p),Q(p))_{p'}\le 4p^{2H/3}, \tag{1}$$
where the subscript removes the full power of $p$. The same statement holds at $p=2$ if each trinomial has three distinct positions. In particular the companion factor $C$ in the seven-unit carry argument satisfies (1) at **every prime**, uniformly in the outer shift. The two independently stripped trinomials in that application cannot be equal: the assigned weight factor $U>1$ would divide the singleton $z(p)$ although $p\nmid U$.

## Two placements suffice

A stripped input supported at only one position gives an immediate bound of three. Otherwise write the two supports as $(0,a,H)$ and $(0,c,d)$ with $0<d\le H$. If $d\le2H/3$, the smaller polynomial already gives the bound. Assume $d>2H/3$. Align the lower endpoints, or shift the second support by $H-d$ to align the upper endpoints. In at least one of these two placements the combined supports have a gap of length at least $H/3$.

Scale to $H=1$. If all gaps in the lower placement were less than $1/3$, reaching $d>2/3$ from zero would require both $a,c$ strictly inside $(0,d)$, with the second of them less than $2/3$. Otherwise some gap in $[0,d]$ is at least $1/3$. Reflect the upper placement across its endpoint to get the same argument with $1-a,d-c$ in place of $a,c$. If both placements failed we would have $1/3<a<2/3$. In the lower placement the first positive position must then be $c<1/3$. In the reflected placement it must be $d-c<1/3$, giving $d<2/3$. Endpoints or positions that coincide cannot evade this argument.

Write $b$ for the upper edge of the selected gap and $g$ for its length. Split the two inputs into $P_0+p^bP_1$ and $Q_0+p^bQ_1$. Multiplication by powers of $p$ preserves the gcd outside $p$. An input on one side alone gives the bound $3p^{H-g}$ outside $p$. Otherwise the common divisor divides
$$ P_0Q_1-Q_0P_1. $$
This is the [determinant at a gap](seven_span_and_gaps.md#the-determinant-at-a-gap), here applied after choosing the better of two separate translations. It is nonzero. If the singleton portions are on the same side, proportionality would make the two original values equal up to a power of $p$. Their base-$p$ digits give identical stripped shapes, since all coefficients are less than $p$ and both constant coefficients are units at $p$. If the singleton portions are on opposite sides, their product is a power of $p$, but a two-unit portion contains a factor greater than one coprime to $p$. Each product has at most four terms and greatest exponent at most $H-g$. Their absolute difference is therefore at most $4p^{H-g}$, proving (1).

## Three tests for a large seven-unit pattern

Let $p=2$ or 3. Consider $N=1+\sum_{j=1}^6p^{a+h_j}$ with $a\ge1$, $0=\min h_j$, $H=\max h_j\ge10^{16}$, and carry-free digit positions. At 2 all six positions are distinct. In a counterexample, if 1 belongs to the triple of the row avoiding $p$, then $p=2$ and the three gaps of the pairs of positive bits are $n,n,n\pm1$ by the [binary alignment theorem](seven_binary_alignment.md). For **every other allocation**, take the binomial $x=X^i(1+X^d)$, trinomial $y=X^j(1+X^e+X^f)$, and singleton $z$ as in the [seven-unit reduction](seven_span_and_gaps.md#which-allocations-need-checking). At least one of the following three conditions must hold.

1. **Two gaps of the trinomial lie on a short list of lines.** For some $(0,t_2,t_3)$ in the [trinomial offset list](binomial_trinomial_gcd.md#the-explicit-factor-in-small-bases), the integer $\gcd(2d,2e-2t_2,2f-2t_3)$ exceeds $H/210$. Thus the three corrected gaps are a common multiple of a vector with absolute entries at most 421.
2. **All four positions of $y+z$ lie on a short list of lines after translation.** Strip their least power and call the remaining three exponents $a_2,a_3,a_4$. For one of the [68 binary or 26 ternary offsets](quadrinomial_small_bases.md), $\gcd(2d,2a_2-\delta_2,2a_3-\delta_3,2a_4-\delta_4)\ge H/(100M_p)$ with $M_2=24$, $M_3=12$.
3. **The binomial gap is three times a common matching gap.** There is a positive integer $w$ with $d=3w$ and $H<49w/11$. The four exponents of $y+z$ split into two pairs whose distances are independently $w$ or $3w$.

These are explicit tests on the six positive exponents. They specialize the [effective codimension-two restriction](effective_sparse_positions.md) to numerical slopes, offsets and matching gaps for seven units. Each lists at most three free integer parameters after the minimum position is set to zero. For the first test these come from one large scale and the starts of three blocks, with one start fixed. For the second there are only two free parameters after translation. The third uses one scale and three block starts before translation, just as in the aligned binary case.

For the proof, the necessary size inequality and (1) give
$$ U^3 A^2\ge (27/16)p^{2H/3}. \tag{2}$$
If $U\ge p^{H/200}$, the [trinomial line proof](binomial_trinomial_gcd.md#bounding-the-remaining-scalar) provides one primitive direction, with a scalar loss at most $p^{6H^{2/3}}$. A nonconstant primitive gcd in this direction requires a listed offset, and has degree at most two. Its roots have modulus $p^{-d/v_0}$ where $v_0$ is its positive binomial slope. It has constant coefficient $\pm1$ since it divides $1+p^d T^{v_0}$ after Laurent monomials are cleared. Its value at 1 is therefore at most $4p^{2d/v_0}$. Primitivity gives $2d/v_0=\gcd(2d,2e-2t_2,2f-2t_3)$. From $1/200-6/H^{1/3}-2/H>1/210$ we obtain the first alternative. These estimates can use the global $H$ in place of the internal block height, with the trivial bound at very small heights.

Otherwise (2) requires $\log_p A>(1/3-3/400)H$. If the four-term estimate has some $W_\delta\ge p^{H/200}$, its definition gives the second alternative. In the remaining case let $k=\deg\gcd(1+X^d,y+z)$. Its polynomial value at $p$ is at most $4p^k$ by the [opposite matchings](seven_large_bases.md#common-roots-at-the-unshifted-pair). Combining this with the [small-base four-term estimate](quadrinomial_small_bases.md) gives
$$ k/H>1/3-3/400-1/200-8/H^{1/4}-2/H>11/35. $$
The odd matching indices start at 3, with $k\le7d/15$. Without index 3 the maximum is $11d/35$, attained at indices 5 and 7. Thus $w=d/3$ is an integer, and both pairs in the index-3 matching have gaps odd multiples of $w$. Since $H<35k/11\le49w/11<5w$, those distances can only be $w$ or $3w$. This proves the third alternative.

Run the exact tests of the two placements, numerical gcds and threshold inequalities from the project directory.

~~~sh
python3 check_two_trinomial_cuts.py
~~~
