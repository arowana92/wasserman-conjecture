# A sharp two-thirds exponent for numerical trinomial gcds

Let $P,Q$ be different, individually normalized sums of three monomials (their lowest exponents are zero). Set $H=\max(\deg P,\deg Q)$, $\epsilon=1/10000$, $\beta=2/3-\epsilon$. For an odd prime $p$ positions may repeat. At $p=2$ assume the positions within each sum are distinct. Then
$$ \gcd(P(p),Q(p))_{p'}\le 4p^{\beta H}\qquad(H\ge10000), \tag{1}$$
with precisely one family exempted from this bound: at $p=2$ the unordered pair of exponent lists
$$ (0,w,3w+2),\qquad (0,2w,3w+1),\qquad w\ge1. \tag{2}$$
Here the notation $p'$ removes the prime power at $p$. For (2) the gcd equals $2^{2w+1}-2^w+1$. Indeed with $Z=2^w$ the two numbers are $(1+2Z)(1-Z+2Z^2)$ and $(1+Z)(1-Z+2Z^2)$, and the linear factors are coprime. Thus no uniform saving in the exponent $2/3$ is possible on this family. For all heights outside (2) one can use $4p^{\beta H+1}$ by the [previous bound](two_trinomial_cuts.md).

## The two shapes near equality

Write the supports as $(0,a,H)$ and $(0,b,D)$ with $D\le H$. A support at a single position gives gcd outside $p$ at most three. Put $T=(1/3+\epsilon)H$. If $D\le H-T$, or the combined supports have a gap of length at least $T$ when either the lower or upper endpoints are aligned, the [two-placement determinant](two_trinomial_cuts.md#two-placements-suffice) gives (1). Assume the opposite throughout this section.

If $D>2T$, the lower placement must have both $a,b$ inside $(0,D)$, and the second of these is less than $2T$. Reflecting the upper placement gives the same conditions on $H-a,D-b$. Thus $H-2T<a<2T$. If both $a,H-a>T$, then $b,D-b<T$, impossible since $D>2T$. Reverse the positions in both supports if needed to take $a\le T$. Now $H-a>T$, so $H-a-T<D-b<T$. In this orientation, with $w=a$,
$$ |H-3w|\le6\epsilon H,\quad |D-b-w|\le3\epsilon H,\quad w\le b\le2T. \tag{3}$$
Indeed $D-b>H-2T$ and $b>D-T>T$. The reversal here also replaces the numerical base by $p^{-1}$, since the original integer evaluations differ from the reversed Laurent evaluations at that base only by powers of $p$.

If $D\le2T$ then $(2/3-\epsilon)H<D\le(2/3+2\epsilon)H$. When $T<a<H-T$ we have $b,D-b<T$. When $a\le T$ the first interior point of the reflected upper placement must be $D-b<T$ and the next is at least $\min(H-a,D)\ge(2/3-\epsilon)H$. Thus $D-b\ge(1/3-2\epsilon)H$ and $b\le(1/3+4\epsilon)H$. Reflecting the argument covers $H-a\le T$. In every case $b$ and $D-b$ lie between $(1/3-5\epsilon)H$ and $(1/3+4\epsilon)H$. With $w=b$ this gives
$$ |D-2w|\le9\epsilon H,\qquad |H-3w|\le15\epsilon H. \tag{4}$$

## Approximating one ratio

Write $B=p$ or $p^{-1}$ according to the orientation. Put $Z=B^w$, $R=B^k$. The two inputs up to powers of $p$ now take the following forms. The first line uses (3) with $k=b$ and the second uses (4) with $k=a$.

| Case | First sum | Second sum | Bounds on offsets |
|---|---|---|---|
| I | $1+Z+B^h Z^3$ | $1+R(1+B^e Z)$ | $|h|\le6\epsilon H,\ |e|\le3\epsilon H$ |
| II | $1+Z+B^e Z^2$ | $1+R+B^h Z^3$ | $|e|\le9\epsilon H,\ |h|\le15\epsilon H$ |

Choose coprime integers $v,u$ with $1\le v\le K=16$, $u\ge0$, $|vk/w-u|\le1/K$ by Dirichlet approximation and division by the common content. Put $s=k-wu/v$. In $\mathbb Q(p^{1/v})[S]$ substitute $Z=S^v$, $R=B^s S^u$, evaluating at $S=B^{w/v}$. Both polynomials have degree at most $4v$, since $u\le3v$ in I and $u\le4v$ in II by the bounds on $k/w$. Let $m$ be their common-factor degree.

For any fixed root $Z$ of the first polynomial in $Z$, the second equation determines $R$ uniquely if it has a solution. In I the vanishing of $1+B^eZ$ would leave the value 1. Since $\gcd(v,u)=1$, the two powers of $S$ determine at most one nonzero $S$. Its multiplicity in the first polynomial is that of $Z$, so $m\le3$ in I and $m\le2$ in II, including repeated quadratic roots.

When $v=1$, in II one has $u\in\{0,1,2,3\}$. A quadratic $1+Z+cZ^2$ with $c>0$ never divides $1+d Z^u+e' Z^3$ with $d,e'>0$. Modulo that quadratic, $Z^2=-(1+Z)/c$ and $Z^3=(1+(1-c)Z)/c^2$. The constant remainder is positive for $u=0,1,3$. For $u=2$ the constant coefficient minus the linear coefficient equals $1+e'/c>0$.

In I, $v=1$ gives $u=1$ or 2. At $u=1$ a quadratic $1+dZ+de'Z^2$ dividing the cubic $1+Z+cZ^3$ would leave a linear quotient with positive constant and leading coefficients, impossible at degree two in the product. At $u=2$, two shared roots counted with multiplicity require a reducible cubic over $\mathbb Q$, since the second cubic has zero linear coefficient and cannot be proportional to the first. Here $c=B^h$ is an integral power of $p$. By the rational-root test a rational root must be $-p^j$ for an integer $j$. The relation $p^j+c p^{3j}=1$ requires its two summands to have the same negative valuation at $p$, and then $p=2$, $p^j=1/2$ and $c=4$. The cubic factors as $(1+2Z)(1-Z+2Z^2)$ with the quadratic irreducible. Using $Z^2=(Z-1)/2$ in $1+dZ^2+de'Z^3$ gives a zero remainder precisely at $d=1$, $e'=2$. These are $B^{k-2w}=1$, $B^{D-k-w}=2$ and $B^{H-3w}=4$. With $B=2$ this gives (2). With $B=1/2$ it would give $D=3w-1>3w-2=H$. Therefore outside (2), $m/v\le3/2$ in I and $m/v\le1$ in II.

## Bounding the numerical gcd

Clear the lowest coefficient power of $p^{1/v}$ separately in each polynomial so its coefficients are algebraic integers. Let $M_1,M_2$ be their degrees and $b_1,b_2$ the spreads of their coefficient exponents of $p$, counting separate summands before collecting equal powers of $S$. Use the [coefficient minor construction](translated_triple_weights.md#a-weighted-determinant-and-an-explicit-bound) over $\mathbb Q(p^{1/v})$. It supplies a nonzero polynomial $J$ of degree $m$ proportional to the gcd, in the integer-ring polynomial span of the two inputs. Each coefficient is a determinant using at most $M_2$ shifts of the first input and $M_1$ of the second. The equal-degree proportional case does not occur since the two supports differ. For any embedding the coefficient absolute values are at most $3^{8v}p^{M_2 b_1+M_1 b_2}$ by the column sums.

The numerical gcd outside $p$ divides $J(B^{w/v})$ in the integer ring with $p$ inverted. If $B=1/p$ multiply this value by $p^{mw/v}$ to make it integral. In both orientations the resulting algebraic integer is nonzero since the two original sums are positive. Every conjugate has magnitude at most $(m+1)3^{8v}p^{mw/v+M_2 b_1+M_1 b_2}$. Taking the norm bounds the numerical gcd by this same quantity, since it is coprime to $p$.

In I the spreads are at most $6\epsilon H$, $3\epsilon H+|s|$, and $M_1=3v$, $M_2\le4v$. In II the corresponding bounds are $9\epsilon H$, $15\epsilon H+|s|$, $M_1=2v$, $M_2\le4v$. Since $v|s|\le w/K$, the bounds on the base-$p$ logarithm become respectively
$$ (3/2+3/K)w+33K\epsilon H+16K+2,\qquad (1+2/K)w+66K\epsilon H+16K+2. $$
Use $w\le(1/3+\epsilon)H$ in I, $w\le(1/3+4\epsilon)H$ in II. Both bounds are less than $\beta H$ for $H\ge10000$. This completes (1).

## Consequence for seven-unit counterexamples

Let $p=2,3$, let $a\ge1$, and suppose $N=1+\sum_{j=1}^6p^{a+h_j}$ has seven carry-free units, with $\min h_j=0$, $H=\max h_j\ge10^{24}$. In a hypothetical counterexample assign $p$ to a nonendpoint row avoiding it. If its units-place 1 belongs to a triple, the [binary alignment theorem](seven_binary_alignment.md) already fixes the three pair gaps to $n,n,n\pm1$. Otherwise let $d$ be the exponent gap of its unshifted pair and $0,e,f$ the relative positions of its unshifted trinomial. Then there are integers $r,s$ with $|r|,|s|\le4$ at 2 and $|r|,|s|\le2$ at 3 for which
$$ \gcd(2d,\,2e-r,\,2f-s)\ \ge\ H/(50000M_p),\qquad M_2=24,\quad M_3=12. \tag{5} $$
Thus all three gaps of these two blocks must agree with integer multiples of one large scale after bounded corrections. The [earlier three tests](two_trinomial_cuts.md#three-tests-for-a-large-seven-unit-pattern) left a third configuration in which the relative start of two matched pairs could move freely. The bound (5) covers that configuration as well.

Use $U,A,C$, the binomial $F$ and four-term sum $y+z$ from those tests. They satisfy $27p^{2H}\le U^3 A^2 C^2$. If $U\ge p^{H/100000}$, the trinomial line bound there gives an allowed offset $(t_2,t_3)$ with $\gcd(2d,2e-2t_2,2f-2t_3)>H/101000$ since $1/100000-6/H^{1/3}-2/H>1/101000$. If some [four-term offset factor](quadrinomial_small_bases.md) has $W_\delta\ge p^{H/100000}$, the definition gives a common divisor $l_\delta\ge H/(50000M_p)$ of $2d$ and the three corrected gaps of $y+z$. Subtracting the corrected exponent of the lowest term of $y$ gives (5). Differences of the doubled heights in the four-term table are at most 4 at 2 and 2 at 3.

Otherwise let $k=\deg\gcd(F,y+z)$, with the monic polynomial gcd evaluated at $p$ bounded by $4p^k$. The four-term estimate gives $\log_p A\le k+\log_p4+H/100000+8H^{3/4}$. The old bound $C\le4p^{2H/3}$ already forces $k/H>11/35$, since $5/200000+8/H^{1/4}+5/H<1/3-11/35$. There is therefore a matching of index 3 and $3\mid d$ by the [matching classification](seven_large_bases.md#common-roots-at-the-unshifted-pair). Also $k\le7d/15$, so $d/H>33/49$.

Case (2) for the two trinomials defining $C$ cannot occur here. For that pair $H\ge3w+2$ and $d$ is a distance between two positions in one of the trinomials. The two largest possible gaps $3w+1,3w+2$ are not divisible by 3. All others are at most $2w+2$. For $w\ge32$ the ratio $(2w+2)/(3w+2)\le33/49$, contrary to the lower bound. For smaller $w$, every gap is below 100 whereas $H\ge10^{24}$.

We may now apply $C\le4p^{\beta H+1}$, allowing the two stripped trinomials to have smaller span than $H$. This gives
$$ k/H>1/3+\epsilon-5/200000-8/H^{1/4}-5/H>1/3+1/20000. $$
A second matching is compulsory. Its odd index $n\ge5$ is coprime to 3 and $d=3ng$, $k=(n+2)g$. The two matchings connect all four positions of $y+z$, with all matched differences divisible by $g$. As $k/d=1/3+2/(3n)>1/3+1/20000$, we have $3n<40000$, $H/g<3(n+2)<40006$. Thus $g\mid d,e,f$ gives (5) with zero corrections. This proves the assertion.

Run the following with python-flint and z3-solver installed. The proof of the real geometric inequalities above is supplemented by an independent exact linear-arithmetic check.

~~~sh
python3 check_sharp_trinomial_gcd.py
~~~
