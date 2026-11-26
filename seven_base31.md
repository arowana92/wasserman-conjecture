# Seven digit units at effective base 31

Let $p$ be prime, $a,b\ge1$, $B=p^b$, and $0=h_1\le\cdots\le h_6=H$. Wasserman's three-row conjecture holds at
$$N=1+p^a\sum_{i=1}^6 B^{h_i}\qquad\text{whenever }B\ge31.$$
Repetitions are allowed at odd primes. At $p=2$ assume the six positions are distinct. In particular **any odd $p\ge31$ dividing $N-1$ with $s_p(N)=7$ suffices**, regardless of its digit positions. This improves the [effective-base cutoff of 64](seven_refined_cutoffs.md).

For odd $p$ the cutoff drops to **$B\ge11$** if some prime $q\mid N$ has $s_q(N)\le5$. The previous sufficient cutoff in this mixed case was 23.

Normalize the smallest position to zero by absorbing a power into the shift. The [span theorem](seven_span_and_gaps.md) allows us to assume $H\ge101$. At odd $p$ a carry among the six positive units leaves at most five units in $N$, which already suffices. Otherwise the [seven-unit reduction](seven_large_bases.md) and [binary alignment](seven_binary_alignment.md) put the 1 in a pair: in the binary case the two required gaps differing by one would both be multiples of $b\ge2$.

Use the notation of the [previous cutoff](seven_refined_cutoffs.md). The two unshifted blocks are $x=X^i F$, $F=1+X^d$, and a trinomial $y$, with singleton $z$ in the shifted pair. The required odd divisors satisfy $U\mid F(B),y(B)$, $A\mid F(B),(y+z)(B)$ and $C\mid y(B),(x+z)(B)$. We have $27B^{2H}\le U^3 A^2 C^2$. Put $h_U=\gcd(F,y)$ and $k=\deg\gcd(F,y+z)$. These polynomial gcds are monic. We use the [signed shift bound](seven_large_bases.md#common-roots-at-the-unshifted-pair):
$$\gcd(F(B),E(B))\le L^{d-j}\gcd(F,E)(B),\qquad j=\deg\gcd(F,E),$$
where $L$ is the Euclidean coefficient norm of $E\bmod F$. In particular $A\le L^{d-k}(B+1)^k$ for $E=y+z$. For odd divisors a common power of two in the reduced coefficients may first be removed.

## Two bounds on the companion factor

**If either $y$ or $x+z$ has three distinct actual positions, $C\le3^H$.** When both do this is the [two-trinomial bound](seven_large_bases.md#bounding-the-factor-from-two-trinomials). Otherwise the repeated-position polynomial, after stripping a monomial, is $X^l+2$ or its reciprocal. The other is $1+X^s+X^f$ with $0<s<f\le H$. A pure monomial times three cannot occur here: such a $y$ either forces the already excluded odd-base scale or has repeated binary digits, and $x+z$ contains the two distinct positions of $x$. Reverse both polynomials if necessary to use $X^l+2$. Resultant absolute values are preserved by reversal. The resultant is nonzero because roots of the monic trinomial with constant 1 are algebraic units, unlike those of $X^l+2$.

For $2l\le H$, Hadamard bounds the resultant by $5^{f/2}3^{l/2}\le(5^{1/2}3^{1/4})^H<3^H$. For $2l>H$, average $|1+X^s+X^f|^2$ over the $l$ roots of $X^l+2$. Cross terms of different residue classes modulo $l$ vanish. Two distinct positions in the same class differ by exactly $l$, so their cross term is negative. The mean square is therefore at most $1+2^{2s/l}+2^{2f/l}\le9^{H/l}$. The geometric mean bounds the resultant by $3^H$.

**If both repeat**, let $l$ be the gap in $y$. The gap in $x+z$ is $d$. By the [two-binomial elimination](seven_refined_cutoffs.md#polynomial-bounds),
$$ C\le 2^{(d+l)/g}+1,\qquad g=\gcd(d,l),\qquad C\le4^H. $$
Also $U\le2^{d/g}+1$: modulo $U$ we have $B^d=-1$ and $B^l=-2^{\pm1}$, so raising to $l/g,d/g$ respectively gives the bound. If $C>3^H$ then $g=1$, since $2^H+1\le3^H$. If $k>0$, coincidence of the two equal terms in $y$ gives $k\mid l,d$ by the [four-term matching argument](seven_refined_cutoffs.md#polynomial-bounds). Thus $C>3^H$ forces $k\le1$.

In fact $C>3^H$ forces $L\le\sqrt6$ for $y+z$ after removing common powers of two. Otherwise the two equal terms of $y$ and $z$ must have the same signed position modulo $2d$. The other coefficient patterns have squared norm at most 6, or are $(2,2)$ up to sign, where two can be removed. Both $d,l>H/2$, since otherwise $2^{(d+l)/g}+1\le2^{3H/2}+1<3^H$ at $H\ge101$. Therefore $z$ actually occurs at the repeated position $r$ of $y$, as $2d>H$. Since $x+z$ also repeats, its repeated position is $r$ as well. The two remaining singleton positions $s,t$ are on the same side of $r$, since their distances $l,d$ exceed $H/2$. The two binomial congruences at $C$ now have the same orientation. Eliminating $B$ gives $C\mid 2^{|d-l|/g}\pm1$. If $d=l$ the two stripped trinomials are identical, already excluded. For different gaps this contradicts $C>3^H$.

## Three distinct residues in the unshifted trinomial

First assume $h_U=1$ and that the three positions of $y$ are distinct modulo $d$. Then $U\le3^{d/2}$ and $C\le3^H$. We have $L\le\sqrt6$ for $y+z$. When $k=0$ this bounds $U^3A^2 C^2$ by $(3^{3/2}\cdot6\cdot9)^H$.

If $0<k\le d/3$ and $L\le2$, the $k$ roots common to $F,y+z$ have $|y|=1$. The mean square of $y$ at the $d$ roots is 3 by orthogonality. Hence
$$ U\le\left((3d-k)/(d-k)\right)^{(d-k)/2}\le4^{(d-k)/2}. $$
Together with $A\le2^{d-k}(B+1)^k$, this bounds $(U^3A^2 C^2)^{1/H}$ by $9\cdot32^{2/3}(B+1)^{2/3}$, since $(B+1)^2>32$.

If $0<k\le d/3$ and $L>2$, the singleton coincides with one of the three signed terms modulo $2d$. At a common root these equal terms must be opposite to both other terms. Interchanging them gives the same matching condition, so the common polynomial is $X^k+1$. At a primitive $2k$-th root all four exponents are therefore congruent modulo $k$. At $k=1$ use the bound $9\cdot3^{3/2}\cdot6(B+1)^{2/101}$ for $(U^3A^2 C^2)^{1/H}$. At $k\ge2$ remove the common monomial from the four terms and use $B^k$ as the base for $U,A$. Their new pair degree is $d/k$, their common degree for $A$ is one, and $B^k+1\le(B+1)^k$. This gives the bound $9(3^{3/2}6)^{1/2}(B+1)^{2/3}$.

It remains that $k>d/3$. All four residues of $y+z$ are distinct. Each nonempty opposite pairing supplies the roots of $X^{d/n}+1$ for odd $n\ge3$. There are two different matchings with distinct indices $n_1<n_2$. Their union gives
$$ k/d=1/n_1+1/n_2-1/\operatorname{lcm}(n_1,n_2). $$
If $n_1=5,n_2=7$, this fraction is $11/35<1/3$. All other pairs with $n_1\ge5$ have $1/n_1+1/n_2\le1/5+1/9<1/3$. Hence $n_1=3$, $n_2=n\ge5$ is prime to 3, and for $g=d/(3n)$ we have $k=d/3+2g\le7d/15$. The two different perfect matchings connect all four positions. Their paired differences are multiples of $d/3$ or $d/n$, so **all four exponents belong to one residue class modulo $g$**.

The monic common factor is $(X^{ng}+1)(X^{3g}+1)/(X^g+1)$ and its value at $B$ is less than $B^k$: for $0<u=B^{-g}\le1/2$ we have $(1+u^n)(1+u^3)<1+u$. Apply the mean-square and signed shift bounds at base $B^g$ to the positions after division by $g$. For $r=k/d$ put $D(r)=4((3-r)/(1-r))^{3/2}$. We obtain
$$U^3A^2\le D(r)^{(d-k)/g}B^{2k}.$$
For $g\ge2$ use $D(r)\le4(19/4)^{3/2}$ and maximize the bound at $r=7/15$. This bounds $(U^3A^2 C^2)^{1/H}$ by $9[4(19/4)^{3/2}]^{4/15} B^{14/15}$. At $g=1$, $k=d/3+2$. If $d\ge300$ then $r\le17/50$, giving the bound $9[4(133/33)^{3/2}]^{33/50} B^{17/25}$. If $d<300$, the verifier uses the exact $r$ at $H=\max(101,d)$, $B=31$, and checks $9^H D(r)^{d-k}B^{2k}<27 B^{2H}$. Its ratio to the right side decreases in both $H$ and $B$.

## Residue collisions and a common polynomial factor

With $h_U=1$ and two occupied residues of $y$, cancellation modulo $F$ would force $U=1$. Otherwise $U\le5^{d/2}$. The two coincident terms agree modulo $2d$. After dividing by a unit we have $B^{s-r}=-2\bmod U$, with the third position $s$ distinct modulo $d$ from the coincident position $r$. Eliminating $B$ using $B^d=-1$ gives $U\le 2^d+1\le3\cdot2^{H-1}$, whether or not $k$ is positive.

For $C\le3^H$ and $k=0$, $L\le\sqrt{10}$ gives the bound $(27/8)^{1/101}\cdot8\cdot10\cdot9$ for $(U^3A^2 C^2)^{1/H}$. For $k>0$ we have $L\le\sqrt6$, $k\le d/3$, and the four positions are all congruent modulo $k$. Indeed, at a primitive root of $X^k+1$ the two coincident terms must be opposite to the two others, as in the preceding section. At $k=1$ the bound is $(27/8)^{1/101}\cdot8\cdot6(B+1)^{2/101}\cdot9$. At $k\ge2$ apply the signed shift estimates at $B^k$, with reduced norms at most $\sqrt5$ and $\sqrt6$ for the two divisors. The latter still holds after removing any common power of two since division of the positions by $k$ preserves signed coincidences. The resulting bound is $9(5^{3/2}6)^{1/2}(B+1)^{2/3}$. If $C>3^H$ then $k\le1$ and $L\le\sqrt6$ by the two-binomial argument. A single bound covers these cases, namely $128(27/8)^{1/101}6(B+1)^{2/101}$. The case when all three terms of $y$ occupy the same signed residue has $U\le3$ and is already [excluded](seven_refined_cutoffs.md#polynomial-bounds) at these bases, including binary with distinct actual positions.

Finally suppose $h_U\ne1$. The [root classification](seven_large_bases.md#the-size-inequality) gives $h_U=X^{2t}-X^t+1$, $d/t=3\bmod6$, $k=0$, $A\le6^{d/2}$ and $C\le3^H$. All three positions of $y$ are congruent modulo $t$. The signed shift bound at base $B^t$ gives
$$U\le3^{(d/t-2)/2}B^{2t}.$$
For $d\ge9t$ and $t\le2$, the resulting bound for $(U^3A^2 C^2)^{1/H}$ is $54\cdot3^{3/2} B^{12/101}$. For $t\ge3$ use $3(d/t-2)/2\le(d-2t)/2$ and maximize at $t/d=1/9$, giving $54\cdot3^{7/18}B^{2/3}$. For $d=3t$ the [earlier argument](seven_refined_cutoffs.md#the-cutoff-at-64) already works at 31: it gives $U<3B^{2t}$ and $A\le3^t+1$, with $C=1$ for $H<6t$ and $(AC)^2\le(9\cdot16^{1/6})^H<B^H$ otherwise.

The 12 uniform bounds from the general cases are strictly less than $B^2$ at $B=31$. Their ratios to $B^2$ decrease thereafter, since each variable factor is a power of $B$ or $B+1$ with total exponent less than two. The verifier raises positive rational powers to a common denominator for the strict comparisons. This proves the general assertion.

## A sparse prime of the weight

Under the mixed hypothesis the [two-base reduction](seven_mixed_primes.md#a-sparse-prime-dividing-the-weight) requires $m=2AC\ge B^H$. For $k=0$ we have $(AC)^{1/H}\le3\sqrt{10}$ if $C\le3^H$, and otherwise $L\le\sqrt6$ gives $4\sqrt6$. For $k=1$ either $C\le3^H$ and $L\le\sqrt6$, or $C>3^H$ and the same norm bound applies. Here $4\sqrt6 (B+1)^{1/101}$ suffices. For $1<k\le d/3$ we have $C\le3^H$. With $L\le2$ the bound is $3\cdot2^{2/3}(B+1)^{1/3}$. For $L>2$ all four exponents are congruent modulo $k$ by the repeated signed term argument. The bound becomes $3\cdot6^{1/4}(B+1)^{1/3}$ after division of the positions by $k$.

When $k>d/3$ the two matchings give $A\le2^{(d-k)/g} B^k$ as above. For $g\ge3$ use $r\le7/15$, giving $(AC)^{1/H}\le3\cdot2^{8/45} B^{7/15}$. At $g=2$ either $n=5,7$ and $d\le42$, or $n\ge11$ and $r\le13/33$. The respective bounds are $3\cdot2^{21/101}B^{98/505}$ and $3\cdot2^{10/33}B^{13/33}$. At $g=1$, $d\ge300$, use $3\cdot2^{33/50}B^{17/50}$. Each of these bounds times $2^{1/101}$ is strictly less than $B$ at 11, and its ratio to $B$ decreases thereafter.

For $g=1$, $d<300$ use the exact index $k=d/3+2$ and check $2\cdot3^H\cdot2^{d-k}B^k<B^H$ starting at $H=\max(101,d)$. The ratio decreases in $B$ and in $H$. At $B=12$ every admissible $d$ passes immediately. At $B=11$ only $H=d=105$ fails. In this case the binomial occupies $0,105$ and one matching pairs the four remaining units at odd multiples of 35. Their residues modulo 105 are distinct, so both gaps are exactly 35. The verifier forms every union of two such pairs in $[0,105]$ and tries each choice of singleton. In all 9,792 cases
$$2\gcd(F(B),(y+z)(B))\gcd(y(B),(x+z)(B))<m,$$
contradicting the required $2AC=m$. This finishes the mixed proof.

Run the verifier with python-flint installed. It checks 20 rational-power inequalities, the small matching indices, 9,792 residual allocations, 7,140 distinct/repeated resultant comparisons, and 48 four-term profiles with common degree exceeding $d/3$ at the specified control degrees.

~~~sh
python3 check_seven_base31.py
~~~
