# Seven units at effective base 64, and a sparse weight prime at base 23

**Theorem.** For an odd prime $p$, positive integers $a,b$, and six nonnegative integers $h_i$ allowing repetitions, the full three-row conjecture holds at
$$N=1+\sum_{i=1}^6 p^{a+b h_i}\qquad (p^b\ge64).$$
In particular every odd predecessor prime $p\ge67$ suffices at digit sum seven, with no restriction on the exponent positions. For $p=2$ the same result holds for six **distinct** $h_i$ and $b\ge6$.

**Mixed theorem.** In the odd-prime assertion the cutoff decreases to $p^b\ge23$ if some prime $q\mid N$ has $s_q(N)\le5$. Here $q$ can be arbitrarily large. In particular for each predecessor prime $p\ge23$ the mixed theorem needs no restriction on the seven digit positions.

## Polynomial bounds

Use the notation of the [seven-unit proof](seven_large_bases.md#the-remaining-allocation). A carry among the six odd-prime powers leaves at most five units. Otherwise the unit at 1 belongs to a pair, by that proof at every odd prime. In binary, the [alignment theorem](seven_binary_alignment.md) forces the same conclusion when $b\ge2$, since two multiples of $b$ cannot differ by one. Normalize $\min h_i=0$, put $B=p^b$, $H=\max h_i$, $m=\sum B^{h_i}$. The inner polynomials are a binomial $x=X^i(1+X^d)$, a trinomial $y$ and a monomial $z$. The assigned scale $U$ and odd factors $A,C$ satisfy
$$U\mid F(B),y(B),\quad A\mid F(B),(y+z)(B),\quad C\mid y(B),(x+z)(B),\qquad F=X^d+1.$$
Here $1\le d\le H$, $U\ge3$, $\gcd(C,B^{3d}-1)=1$, and $27B^{2H}\le27m^2\le U^3 A^2 C^2$. The first proof gives $h=\gcd(F,y+z)$ monic of degree $k\le7d/15$, and its signed shift matrix bound
$$ A \le L^{d-k}(B+1)^k, $$
where $L$ is the Euclidean norm of the coefficients of $y+z$ modulo $F$. Both $A,C$ are odd, so a factor 2 in a reduced coefficient vector may be removed before applying the bound.

**The factor $C$.** In all cases $C\le4^H$. If either stripped trinomial $y,x+z$ has three distinct exponents, the [Sylvester bound](seven_large_bases.md#bounding-the-factor-from-two-trinomials) gives $C\le15^{H/2}<4^H$ (and $3^H$ if both do). Otherwise they are binomials with coefficient magnitudes 1 and 2, at respective distances $l,d$. With $g=\gcd(l,d)$, eliminating $B$ from the two binomial congruences makes $C$ divide a nonzero integer $2^v\pm1$ with $v\le(l+d)/g$. Indeed the equations are $B^l=-2^{\pm1}$ and $B^d=-2^{\pm1}$ modulo the odd $C$. The zero identity requires $l=d$ and the same orientation, in which case the stripped polynomials agree, already impossible since $U$ would divide $z(B)$. For $l\ne d$, $v\le2H-1$. For $l=d$ with opposite orientations $C\mid3$.

If the three terms of $y$ are in distinct residue classes modulo $d$, then **$C\le3^H$** even when $x+z$ has a repeated position. In the latter case the stripped $x+z$ is $X^d+2$ or its reciprocal. Reverse both stripped polynomials for the reciprocal case, which preserves common divisibility after evaluation at the unit $B^{-1}$ modulo $C$. Let the other stripped polynomial have exponents $0,s,l$ with $l\le H$. At the roots of $X^d+2$ orthogonality gives mean square $1+2^{2s/d}+2^{2l/d}\le1+2\cdot4^{H/d}\le9^{H/d}$; the final bound uses $(1+4+4)^t\ge1+4^t+4^t$ for $t\ge1$. The nonzero resultant is bounded by $3^H$ by the geometric mean. The resultant is nonzero because the roots of the monic positive trinomial with constant coefficient one are algebraic units, and those of $X^d+2$ are not.

**Residue collisions.** If two of the four positions in $y+z$ coincide modulo $d$, then $k\le d/3$. For agreement modulo $2d$, opposite pairing at a root forces the two equal terms to match the other two, giving only one root condition after interchanging equal terms. If the terms differ by $d$ modulo $2d$, they cancel at every root of $F$, and the other two must cancel. In either case a nonempty proper root set is the root set of a single $X^k+1$ with $d/k$ odd and at least three. This includes both signs in a reduced binomial. In the cancellation case the differences used for the opposite pairing are the original exponent differences, even if cancellation leaves a signed binomial.

If two terms of $y$ agree modulo $2d$ at position $r$ while the third has a different residue modulo $d$ at $s$, a nonempty such root set forces $k\mid s-r$. At every common root the third term must be opposite to one of the two equal terms; the matching description gives the divisibility. As $2+B^{s-r}=0\bmod U$, raising to $d/k$ gives
$$U\le2^{d/k}+1\quad(k>0).$$
In this case if $k>0$ the reduced vector for $y+z$ has norm at most $\sqrt6$ after removing any common factor 2. In fact the three occupied entries can have magnitudes $2,1,1$. With only two entries the possible larger patterns are $2,2$ (remove two) and $3,1$ (no unit-circle zero). Opposite signs can only decrease the magnitudes. If two terms of $y$ cancel modulo $F$, then $U=1$. If all three have the same signed residue, $U\le3$, $A\le4^d$ and $k=0$. The latter case is impossible at an odd base, where $U\ge7$, and gives $U^3A^2C^2\le27\cdot144^H$ in binary, since all six positions there are distinct.

## The cutoff at 64

First assume $\gcd(F,y)=1$. When $y$ has three distinct residues modulo $d$, $U\le3^{d/2}$. If $k\le d/3$ use $L\le\sqrt6$. For $k>d/3$ all four residues of $y+z$ are distinct by the collision argument, so $L=2$. The resultant of $F,y$ is nonzero, and at each of the $k$ roots of $h$, $y=-z$ has modulus one. The sum of $|y|^2$ at all $d$ roots of $F$ is $3d$. Thus the geometric mean on the other roots gives
$$U\le ((3d-k)/(d-k))^{(d-k)/2}\le(19/4)^{(d-k)/2}.$$

For a collision between two terms of $y$ use $U\le5^{d/2}$ and $L\le\sqrt{10}$ when $k\le d/9$. For $k>d/9$ we have $U\le129$ since $d/k=3,5,7$, and $k\le d/3$. In the following table each entry bounds $(U^3A^2 C^2/B^{2H})^{1/H}$ from above. In the last row take $H\ge21$. The estimates follow by using $d\le H$ and maximizing the powers of $(B+1)$ at the stated upper bounds on $k/d$.

| Case | Upper bound |
|---|---|
| Distinct residues of $y$, $k\le d/3$ | $3^{3/2}\,9\,6^{2/3}(B+1)^{2/3}/B^2$ |
| Distinct residues of $y$, $k>d/3$ | $(19/4)^{4/5}\,9\,4^{8/15}(B+1)^{14/15}/B^2$ |
| Two residues of $y$, $k\le d/9$ | $5^{3/2}\,16\,10^{8/9}(B+1)^{2/9}/B^2$ |
| Two residues of $y$, $k>d/9$ | $129^{3/21}\,16\,10^{2/3}(B+1)^{2/3}/B^2$ |

All four expressions decrease with $B$ and are less than one at 64 by integer powers. For $H\le20$ the [21-position theorem](seven_progressions.md) already covers odd bases. For binary, or alternatively for both, the last collision case at these heights satisfies
$$\min(5^{3d/2},(2^{d/k}+1)^3)\,10^{d-k}(B+1)^{2k}\,16^H < 27 B^{2H}.$$
The verifier checks every $1\le d\le H\le20$ and $1\le k\le d/3$ with $d/k$ odd integral, at $B=64$. The quotient by $B^{2H}$ is decreasing. This proves all the coprime cases.

For $\gcd(F,y)\ne1$, use the [root classification](seven_large_bases.md#the-size-inequality). The common factor is $X^{2t}-X^t+1$ with $d/t=3\bmod6$, the three residues of $y$ are distinct modulo $d$, and $k=0$. If $d\ge9t$, the signed shift bound gives $U\le3^{(d-2t)/2} B^{2t}$. Also $A\le6^{d/2}$ and $C\le3^H$. The ratio in the table is now at most
$$3^{3/2}\cdot54\,(B^6/27)^{1/9}/B^2<1\qquad(B\ge64).$$

If $d=3t$, the earlier proof gives $U<3B^{2t}$ and $A\le3^t+1$. When $H<6t$, the three distinct exponents of $y$ occupy one class modulo $2t$ and have span less than $6t$, hence are exactly $j,j+2t,j+4t$. The stripped trinomial $1+X^{2t}+X^{4t}$ divides $X^{6t}-1$. Since $C$ divides $y(B)$ and is coprime to $B^{3d}-1$ and to $F(B)$ (as $C\mid(x+z)(B)$), it is coprime to $B^{6t}-1=(B^{3t}-1)F(B)$. Hence $C=1$. As $H\ge4t$ and $3^t+1<B^t$, the required size fails. When $H\ge6t$ use $C\le3^H$ and $(3^t+1)^2\le16^t$. This gives $(AC)^2\le(9\cdot16^{1/6})^H<B^H\le B^{2H-6t}$, again a contradiction. This completes the first theorem.

## The cutoff at 23 with a sparse weight prime

By the [mixed obstruction](seven_mixed_primes.md#a-sparse-prime-dividing-the-weight), only $N=3\bmod4$ and $s_q(N)=5$ remain, and the companion row must carry at every odd predecessor prime. Thus every odd prime power in $m$ is assigned to the $p$-avoiding row and $m=2AC$. We may assume $H\ge21$ by the 21-position theorem. We need only bound $(2AC/B^H)^{1/H}$, independently of the gcd of $F$ with $y$.

For distinct residues of $y$ use $C\le3^H$, and $L\le\sqrt6$ if $k\le d/3$ or $L=2$ if $k>d/3$. The respective bounds are
$$2^{1/21}\,3\,6^{1/3}(B+1)^{1/3}/B,\qquad
2^{1/21}\,3\,2^{8/15}(B+1)^{7/15}/B.$$
For two residue classes use $C\le4^H$. If $k=0$, $L\le\sqrt{10}$ gives the bound $2^{1/21}4\sqrt{10}/B$. If $k>0$, $L\le\sqrt6$ after removing any factor 2, and $k\le d/3$. This gives $2^{1/21}4\,6^{1/3}(B+1)^{1/3}/B$. Each expression decreases with $B$ and is less than one already at 23. This contradicts $2AC=m\ge B^H$.

At the remaining primes there is also a cofactor test. Put $Q=q^c$ for the full weight power and let $D$ be the highest exponent in the base-$q$ expansion of $N/Q$. A counterexample requires **$1\le D$ and $c\ge2D+1$**. In particular the mixed theorem holds for every predecessor prime if $q^3\nmid N$.

Indeed the companion row avoiding $q$ has a pure $q$-power part. Its assigned weight scale can only be $Q$, and it must take the full power since the other row cannot avoid $q$. The $p$-avoiding row carries at 2 and must take every odd predecessor power, with product $(N-1)/2$. Its normalized parts sum to $Q$ and their product is divisible by $((N-1)/2)^2$. Thus $Q^3\ge27(N-1)^2/4\ge3N^2$, so $Q\ge3(N/Q)^2>q^{2D}$. The two even parts of the companion have two $q$-units each and opposite gap parities; hence the larger occupied exponent strictly exceeds $c$. This proves the test, including the case of a cubefree weight.

## Verification

Run `python3 check_seven_refined_cutoffs.py` with python-flint. It checks the numerical bounds as rational integer powers, including 163 small-height triples, and compares 28,140 exact Sylvester determinants with the bounds from the binomial shapes. Four-position checks include 6,428 proper collision factors, 6,600 two-residue scale comparisons, and 20 cases with common degree exceeding $d/3$ and distinct residues. Direct enumerations at the improved bases check 83,046 carry-free rows. A separate comparison at heights 21, 30, 42 and 63 tests 402,492 rows at 2,812 mixed weight/base/shift choices, including 719 row scales of at least seven. [Verifier](check_seven_refined_cutoffs.py)
