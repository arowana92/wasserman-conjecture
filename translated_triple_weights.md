# Finitely many inner exponent triples for a translated three-term sum

Let $p$ be prime, $a,k\ge1$, and $0<i<j$. Put
$$N=1+p^a(1+p^i+p^j)(1+p^k),\qquad H=j+k.$$
**The full three-row conjecture holds if $H>414{,}720$.** At $p\ge5$ the upper bound for an exceptional $H$ decreases to $262{,}656$. Thus only finitely many quadruples $(p,i,j,k)$ can admit a counterexample at any outer shift. The finite interval between the known small-span theorem and the upper bounds has not been exhaustively checked.

**Bounded column span.** For $0<i<j\le48$ the conjecture holds for every prime and all positive integers $a,b,k$ at
$$1+p^a(1+B^i+B^j)(1+B^k),\qquad B=p^b.$$
The translation $k$, outer shift $a$ and dilation $b$ are unrestricted. The [arithmetic-progression theorem](geometric_product_weights.md) handles $j=2i$ without a bound on $j$.

**Independent column and layer spacings.** The conjecture also holds for every $0<s<t\le5$ at
$$1+p^a(1+p^{sd}+p^{td})(1+p^k)$$
with $a,d,k$ arbitrary positive integers. Here the layer translation need not be divisible by the column spacing $d$.

## Binary mergers and the location of the unit at 1

Write $m=(1+p^i+p^j)(1+p^k)$. At odd $p$ each position occurs at most twice. The [seven-unit reduction](seven_large_bases.md) puts the units-place 1 in a pair in the row avoiding $p$ in any counterexample. At $p=2$ suppose first the six positions in $m$ are distinct and 1 belongs to a triple. [Binary alignment](seven_binary_alignment.md) requires three disjoint pairs of the six positions of gaps $n,n,n\pm1$, with $27m^2\le(2^n+1)^3$. Since $m>2^H$ and $2^n+1\le 3\cdot2^{n-1}$, we have $2H+3<3n$. All three gaps exceed $H/2$, so each pair crosses from the bottom three positions to the top three. If $k>j$, their lengths sum to $3k$, impossible for $3n\pm1$. If $k<j$, the bottom three are $0,i,k$ and the top three are $j,i+k,j+k$. The sum is $2j+k$. This gives $3n\le2j+k+1\le2H$, also impossible.

A binary repetition requires $k=i,j$, or $j-i$. Merging one pair of equal powers gives these respective lists of positions:
$$
(0,i+1,2i,j,j+i),\quad (0,i,j+1,i+j,2j),\quad
(0,i,j-i,j+1,2j-i).
$$
A further merger leaves at most five units including the outer 1, already [covered](six_digit_gaps.md). Otherwise the [six-unit equal-gap theorem](equal_predecessor_gaps.md) requires four distinct entries forming two pairs of the same gap $n$, with $27m^2\le(2^n+1)^3$ as before. For those three lists the integral inequalities $0<i<j$, $2H+3<3n$ force $j=2i$, $k=j$, $n=3i$. The supplied [Lean proof](TranslatedTripleGaps.lean) exhausts the omitted position, pairing and gap orientations using linear integer arithmetic. The [proof for $j=2i$](geometric_product_weights.md#coincident-positions) excludes this last case. Thus throughout the remaining proof the units-place 1 is paired with a singleton from $m$, and the six positions in binary are distinct.

## Normalizing by the singleton

The row is $(p^a x_0,p^a y_0,1+p^a z_0)$ with $x_0,y_0,z_0$ sums of two, three and one labelled terms of $m$. The [prime-power assignment](seven_span_and_gaps.md#which-allocations-need-checking) requires odd $U,A,C$, coprime to $p$, with
$$ U\ge3,\quad (U,m)=1,\quad U\mid x_0,y_0,\quad
 A\mid x_0,y_0+z_0,\quad C\mid y_0,x_0+z_0,\quad 27m^2\le U^3(AC)^2. $$
In particular $S=AC$ divides both $m$ and $x_0y_0$ since $A,C$ are coprime.

Divide the six monomials by $z_0$. Write the resulting grid as $1,X,Y,Z,ZX,ZY$, where $1,X,Y$ are in the singleton's layer, and $Z=p^{\pm k}$. Its exponent vector $E=(\log_p X,\log_p Y,\log_p Z)$ has height
$$ h(E)=\max(0,E_1,E_2)-\min(0,E_1,E_2)+|E_3|=H. $$
In these signed coordinates put $x=x_0/z_0$, $y=y_0/z_0$, $M=m/z_0=(1+X+Y)(1+Z)=1+x+y$. Divisibilities away from $p$ hold in $\mathbb Z[1/p]$. Up to exchanging $X,Y$, only three choices for $x,y$ need checking.

| Type | $x$ | $y$ |
|---|---|---|
| I | $X+ZY$ | $ZX+Y+Z$ |
| II | $X+Z$ | $ZX+Y+ZY$ |
| III | $Z(1+X)$ | $X+Y+ZY$ |

A pair in the same column divides $M$ up to a monomial. The pair $X+Y$ leaves $y=Z(1+X+Y)$, also a divisor up to a monomial. Both force $U\mid m$. The pair $Z(X+Y)$ leaves $y=X+Y+Z$, forcing $U\mid Z$. This exhausts the ten ways to choose two of the remaining five grid entries.

## Common roots on an affine line of exponents

On an integer line $E=w+n v$ with $v$ primitive, substitute $(X,Y,Z)=(p^{w_1}T^{v_1},p^{w_2}T^{v_2},p^{w_3}T^{v_3})$. Write $h=h(v)$. From each polynomial separately strip its lowest power of $T$ and its lowest coefficient power of $p$ among its monomials. Denote the resulting integer polynomials by the same letters. Common factors below are taken over $\mathbb Q$ with degrees
$$d_u=\deg\gcd(x,y),\qquad d_s=\deg\gcd(M,xy).$$
The grid gives $\deg M=h$ and all stripped polynomials have positive constant term. All identities between the original expressions can be used on the complex torus $T\ne0$.

We claim that on every **nonhomogeneous** line containing a feasible exponent vector,
$$2h-3d_u-2d_s\ \ge\ h/3. \tag{1}$$
Here feasibility only uses the distinctions between the three columns, $Z\ne1$, the six distinct positions at $p=2$, and $U>1$ dividing $x,y$. At $p\ge5$ the lower bound is $h$.

To see this, consider three systems of common zeros:
$$ (x,y),\qquad (x,1+X+Y),\qquad (y,x+1). $$
Call the associated degrees $d_u,d_l,d_c$. The remaining factor of $M$ contributes at most $h/2$ to $\deg\gcd(x,M)$. In fact at $Z=-1$ the binomial $x$ becomes a signed sum of two distinct column monomials. The degree of its common factor with $1+Z$ is at most the smaller of $|v_3|$ and the difference in slopes of those columns in absolute value. If $v_3=0$, $1+Z$ is a nonzero constant. If the column slopes agree, cancellation would require identical column magnitudes on the entire line, contrary to distinctness. For nonzero slopes the bound follows from the two binomial equations, which have simple roots. Consequently
$$d_s\le h/2+d_l+d_c,$$
using $\gcd(y,M)=\gcd(y,x+1)$ on the torus and bounding multiplicities in the product $xy$.

The table gives a character value and a vanishing three-term sum at each possible common root. A character here is just a monomial with integer exponents.

| System | I | II | III |
|---|---|---|---|
| $x=y=0$ | $X/(ZY)=-1$, $y=0$ | $X/Z=-1$, $y=0$ | $X=-1$, $y=0$ |
| $x=1+X+Y=0$ | $X/(ZY)=-1$, $1+X+Y=0$ | $X/Z=-1$, $1+X+Y=0$ | empty |
| $y=x+1=0$ | $Z=\pm1$, $1+X+ZY=0$ | $Y/Z=1$, $1+X+Z=0$ | $YZ=1$, $y=0$ |

For the bottom row take $y-Z(x+1)$ in types I and II and $Zy-(x+1)$ in III. The respective expressions are $Y(1-Z^2),(1+Z)(Y-Z),(1+Z)(ZY-1)$. In II and III, $Z=-1$ would force $X=0$. No line can satisfy a listed character equation identically and contain a feasible vector. A negative value is impossible at positive inputs. The value $Z=1$ is excluded by $k>0$. In II the identity $Y=Z$ would give $y=Z(x+1)$, and in III $YZ=1$ would give $Zy=x+1$, both incompatible with $U>1$. Thus each common root is simple in the gcd, by the corresponding binomial identity. In the bottom row the factors removed at a common root are local units.

Every root $r$ satisfies $r^d=\pm p^l$ with a nonzero integer $d$. In particular $r=\zeta p^\lambda$ for a root of unity $\zeta$ and a rational $\lambda$. For any one character equation its roots have a single absolute value. At that absolute value, the three complex summands in its listed trinomial have fixed magnitudes. Their relative phases have at most two possibilities: a triangle with prescribed side lengths is fixed up to rotation and conjugation, including the collinear case. In the phase coordinates of $X,Y,Z$, the character and the two ratios to the first trinomial term define a monomial map of determinant of absolute value 1 except for the top left cell, where it is 2. These determinants follow by subtracting the three-component exponent vectors in the table. Each admissible phase triple has at most that many preimages by integer diagonalization. The map $r\mapsto (p^{w_i}r^{v_i})_{i=1}^3$ is injective since $v$ is primitive. Accounting for the two signs in the bottom left cell gives
$$ (d_u,d_l,d_c)\le (4,2,4)\ \text{in I},\qquad (d_u,d_l,d_c)\le(2,2,2)\ \text{in II and III}. $$
These bounds prove (1) for $h\ge36$ in I or $h\ge21$ in II and III.

To treat small directions and $p\ge5$, put $t=w+\lambda v$. A summand at the root is a root of unity times $p^s$ with $s$ the scalar product of its exponent vector and $t$. At a place over $p$ the three valuations are these same $s$, normalized by $v_p(p)=1$. Thus the two least $s$ agree. If the third is larger by $\gamma>0$, squaring the absolute value relation gives
$$p^{2\gamma}=2+\omega+\omega^{-1}$$
for a root of unity $\omega$. Every conjugate of the right hand side is real in $[0,4]$. Writing $2\gamma=m_1/n_1$ in lowest terms, $\mathbb Q(p^{m_1/n_1})=\mathbb Q(p^{1/n_1})$ by Bézout, of degree $n_1$ by Eisenstein. For $n_1>1$ the conjugates of $p^{m_1/n_1}$ include a negative or nonreal value. Therefore $p^{2\gamma}$ is one of $2,3,4$. Equal magnitudes in the trinomial force $t=0$ by the nonsingular character matrix, putting the line through the origin. On a nonhomogeneous line only $p=2$ with $\gamma=1/2,1$, or $p=3$ with $\gamma=1/2$, can occur. In particular at $p\ge5$ all three degrees vanish and $2h-3d_u-2d_s\ge h$.

For $p=2,3$ the [exact degree verifier](check_translated_triple_lines.py) enumerates $t$ by putting the character height equal to zero and assigning the three summand heights, up to a common additive constant, as a permutation of $(0,0,\gamma)$. All coordinates of $4t$ are integral. For every primitive direction in range, a line through $t$ has an integer representative $t+s v/4$ for some $s=0,1,2,3$ if it has any integer point, by primitivity. After canonicalizing the origin the verifier independently computes $d_u,d_s$ from the polynomials. It tests 745,332 oriented lines. Exactly two checks fall below (1), the two orientations of one binary line on which $ZX=Y$. The binary merger proof covers it. On lines containing none of the listed offsets there is no common root of the three systems, so $d_s\le h/2$ suffices. This proves (1).

The independent [C++ line verifier](check_translated_triple_lines_cross.cpp) generates the possible offsets by checking the triangle heights on a quarter-integer grid. The offsets from the exact linear systems range from $-2$ to $2$ in each coordinate. The verifier keeps one orientation and computes degrees after reduction modulo 1009. The leading coefficients stay nonzero in each accepted computation, so both modular gcd degrees bound their rational counterparts from above. All 361,811 eligible incidences pass at this first modulus.

## A weighted determinant and an explicit bound

For polynomials $f,g$ of degrees $l,m'$ with common-factor degree $d$, use the $s=l+m'-2d$ columns $T^i f$ ($0\le i<m'-d$) and $T^i g$ ($0\le i<l-d$). When $s>0$, they are independent by removing the gcd and using the degree bounds. Over $\mathbb Q$ they span all multiples of the gcd of degree at most $s+d-1$, by dimension. The $s-1$ coefficient rows above degree $d$ have rank $s-1$. Their signed maximal minors give a combination of the columns that is nonzero and has degree $d$. Call it $R$. Each coefficient of $R$ is, up to sign, the determinant obtained by adjoining its coefficient row to those $s-1$ rows. If $s=0$, use $R=f$. In either case $R$ is a scalar multiple of the gcd and is nonzero at $T=1$ for our positive-coefficient pairs.

Use this construction on the line through an actual exponent point $E$, with direction $v$ primitive. Substitute $p^{E_i}T^{v_i}$ and strip only the lowest monomial power of $T$ from each polynomial. These operations take place in $\mathbb Z[1/p][T]$, so each assigned divisor divides $R(1)$ away from $p$. Write $E=\alpha v+r$ for a **real** $\alpha\ge0$, and set $D=h(r)$, $h=h(v)$. Coefficient bounds follow by keeping track of individual grid monomials in the determinants. For a term from row $t$ of a column $T^i f$, its $v$-weight is $t-i$ plus the lowest original weight in $f$. Therefore all monomials of the determinant for coefficient $j$ have the same weight $W+j$ with $W$ independent of $j$. The entire spread of the $\alpha v$ contribution to the powers of $p$ in $R(1)$ is at most $\alpha d$.

The residual $r$-weights of the six grid monomials have spread $D$. For the $x,y$ pair a determinant has at most $h$ columns of each kind, giving residual spread at most $2hD$. The sum of absolute multiplicities in one coefficient is at most $2^h3^h$: expand the determinant and bound its unsigned sum by the product of the column sums of multiplicities. For $M,xy$ the column counts are at most $2h,h$, with six monomials per column in either case and residual spreads at most $D,2D$. This gives multiplicity at most $6^{3h}$ and spread at most $4hD$. The $s=0$ choices satisfy the same bounds since $h\ge1$. There are at most $h+1$ coefficients in each $R$. Removing the smallest power of $p$ from its nonzero value leaves an integer divisible by the assigned divisor. Consequently
$$ U\le(h+1)6^h p^{\alpha d_u+2hD},\qquad S\le(h+1)6^{3h}p^{\alpha d_s+4hD}. \tag{2} $$
This estimate uses real approximation only in the size argument. The determinants and divisibility remain integral away from $p$.

Since $H\ge\alpha h-D$, $m\ge p^H$ and $27m^2\le U^3S^2$, taking logarithms in (2) gives
$$ \alpha(2h-3d_u-2d_s)\le h(30+16D). \tag{3} $$
Here $h+1\le p^h$, $6^h\le p^{3h}$, $6^{3h}\le p^{8h}$. For $p\ge5$ the constant 30 decreases to 19 using $6^h\le p^{2h}$ and $6^{3h}\le p^{4h}$.

Choose the approximation in the original coordinates $e=(i,j,k)$ and put $K=\max(j,k)$. For an integer $Q\ge2$, there are $1\le q\le Q^2$ and an integer vector $v_0$ of norm $q$ satisfying $\|qe/K-v_0\|_\infty\le1/Q$. One coordinate of $e/K$ is 1. Put the fractional parts of $0,1,\ldots,Q^2$ times the other two coordinates in $Q^2$ square boxes and subtract two in the same box. The maximal coordinate is represented exactly. Write $v_0=g v'$ with $v'$ primitive and set $\alpha=gK/q$. The column spread plus the absolute layer difference of $e-\alpha v'$ is at most $2K/(qQ)$. If the exact coordinate is $j$ the column spread is at most $K/(qQ)$. If it is $k$ the layer difference is zero. Relabelling columns relative to the singleton and reversing layers if needed is an integral unimodular change of coordinates, preserving this height for real vectors. It takes $v'$ to a primitive $v$ in the coordinates of (2). Thus $D\le2\alpha/Q$.

On a homogeneous line $e=bv'$ with $b$ integral. Then $K\le5q$ or $|b|\ge6$. In the latter case all relative positions are multiples of $|b|$ and the [effective-base theorem](seven_refined_cutoffs.md) applies since $p^{|b|}\ge64$. For nonhomogeneous lines take $Q=48$. When the first three common-root systems in the table are empty on this line, $d_u=d_l=d_c=0$, giving a margin of at least $h$ even in the small bases. Hence $\alpha\le30+16D\le30+32\alpha/Q$, so $\alpha\le90$. This proves $H\le2K\le2\alpha Q^2\le414{,}720$ in this case. At $p\ge5$ these systems are always empty on a nonhomogeneous line and the constant 19 gives $H\le2\cdot57\cdot48^2=262{,}656$.

If one system has a common root, the nonhomogeneous line passes through one of the listed offsets $t$. Their heights are at most 3, as verified from the same character matrices. Write $E=t+\beta v$. In original coordinates the component attaining $K$ has positive slope and the corresponding offset coordinate has absolute value at most 3, since column differences and the layer spacing are each bounded by the offset height. Thus $\beta>0$ for $K>3$. Apply (3) with this exact offset to get $\beta(2h-3d_u-2d_s)\le78h$. If $h\le48$, (1) gives $H\le234h+3\le11{,}235$. Otherwise the degree bounds from the table give margin at least $h-24$. In the original coordinates all entries of $v_0$ are in $[0,q]$, so $h\le2q\le4608$. Since $h^2/(h-24)\le h+48$ for $h\ge48$, we obtain $H\le78(h+48)+3\le363{,}171$. This proves the uniform bound in all cases. Finally all primes $p\ge67$ are already covered at every span by the effective-base theorem.

## Counting the remaining outer shifts

There are effective constants $C$ and $0<\theta<1$ such that in **any** interval of $L\ge2$ consecutive outer exponents $a$, at most $C L^\theta$ choices of $(p,i,j,k,a)$ from the entire family can fail. Only the finitely many inner quadruples within the height and prime bounds need be considered.

Fix one allocation at such a quadruple and retain the unscaled positive integers $x_0,y_0$. The assigned $U$ has an odd prime divisor $q\mid\gcd(x_0,y_0)$ with $q\ne p$, and this row must avoid $q$. Therefore the sum of $p^a x_0$ and $p^a y_0$ has no base-$q$ carry. Choose $b$ with $q^{b-1}\ge\max(x_0,y_0)$. If the input $p^a$ has $b$ consecutive digits $q-1$, the last digit of that block in either product is $q-1$. Indeed for multiplier $d$ the incoming multiplication carry $c$ satisfies $0\le c<d$. The output block modulo $q^b$ is $c-d\bmod q^b$, between $q^b-q^{b-1}$ and $q^b-1$. The two products would carry.

Take $l$ such that $q^l\le L<q^{l+1}$, allowing $l=0$. Among the $q^l$ residue strings of length $l$, at most $q^b(q^b-1)^{\lfloor l/b\rfloor}$ avoid the forbidden block, by testing just the disjoint aligned blocks. Each residue occurs $O_{p,q}(1)$ times in an interval of $L$ exponents. In fact, for $d_0=\operatorname{ord}_q(p)$ and $c_0=v_q(p^{d_0}-1)$, the order modulo $q^l$ is $d_0 q^{l-c_0}$ for $l\ge c_0$, by odd-prime lifting. Smaller $l$ are bounded. This gives the exponent $\log(q^b-1)/\log(q^b)<1$ for the allocation and prime $q$. Taking a finite sum and the largest exponent proves the claim.

## The bounded column span calculation

The [span-100 theorem](seven_span_and_gaps.md) covers normalized $j+k\le100$ with arbitrary dilation. For $k>j$ its gap bound at the gap $k-j$ requires $B^{2k-12j}\le41472$, so $k\le6j+7$ in any exception. Binary coincidences at any dilation have been handled above. It suffices to test the 26 prime powers $B<64$, $1\le i<j\le48$, and $\max(1,101-j)\le k\le6j+7$.

In each of the 36 remaining labelled allocations the two blocks and the singleton are linear in $T=B^k$. For two blocks $f=f_0+f_1T,g=g_0+g_1T$ a nonzero determinant $f_0g_1-f_1g_0$ bounds their gcd away from $p$. A one-term block supplies its coefficient after removing a unit power of $T$. Powers of $2$ and $p$ can be removed from all caps. The verifier first applies these caps to the three pairs for $U,A,C$, and directly tests the remaining integers using $U$ bounded by the part of $\gcd(x_0,y_0)$ coprime to $2pm$ and $AC$ bounded by the odd part of $\gcd(m,x_0y_0)$.

At normalized width 48, 2,148,696 numerical size tests remain with a sufficiently large scale. None requires a companion search. An earlier run preserving all 60 allocations and retaining the base factors in the determinant caps also excluded every case.

## Independent column and layer spacings through five

Reduce $s,t$ by their gcd, absorbing it into $d$. The ratio $1:2$ is [proved separately](geometric_product_weights.md). For the eight reduced ratios through $t=5$ the [sector verifier](check_translated_triple_ratios.py) uses the same unimodular subdivisions as the [geometric-sum proof](geometric_product_weights.md#resultants-on-a-sector), with column coordinates $0,s,t$. Its linear height for the square of $m$ is $L(v)=2(tv_1+v_2)$. In each oriented sector write $(d,k)=(b+\ell)v+\ell w$. For the three assigned divisor pairs project away $T=p^{b+\ell}$, leaving integer polynomials $D_j(Y)$ in $Y=p^\ell$. Each is nonzero and bounds its divisor by its absolute evaluated value when that value is nonzero. A block independent of $T$ supplies the polynomial directly. As before, strip monomials and otherwise use the integral resultant. The recorded subdivision verifies
$$g=L(v+w)-\sum_{j=0}^2 u_j\deg D_j>0,\qquad (u_0,u_1,u_2)=(3,2,2).$$
Away from the roots, a counterexample needs $27p^{bL(v)+\ell g}\le\prod_j\|D_j\|_1^{u_j}$, which gives explicit finite loops.

At a root $Y=1$, $\ell=0$, so $(d,k)=bv$. Here $p^b<64$ suffices by the large-base and binary merger theorems. Every other positive integer root fixes $p,\ell$ if it is a prime power. Specialize along this line and compute $d_u,d_s$ directly from $x,y,M$. The direction is primitive in the three exponent coordinates: its original vector is $(sv_1,tv_1,v_2)$, with both $(s,t)$ and $(v_1,v_2)$ coprime. Its height is $h=tv_1+v_2$ and the residual height is $D=\ell(tw_1+w_2)$. The program verifies $2h-3d_u-2d_s>0$ in each such case, so (3) explicitly bounds $b+\ell$. Across 17,352 oriented sectors the only positive integer roots are $1,2,3,8$. There are 156 nonhomogeneous root checks and 230,276 direct integer tests after pruning large effective bases and binary repetitions. Every test violates the necessary size inequality.

## Reproduction

Use Python with SymPy and python-flint and Lean 4.12.0. Run these from the repository root.

~~~sh
lean TranslatedTripleGaps.lean
python3 check_translated_triple_lines.py
c++ -O3 -std=c++17 check_translated_triple_lines_cross.cpp -o /tmp/translated-triple-check
/tmp/translated-triple-check
python3 check_translated_triple_columns.py
python3 check_translated_triple_ratios.py
~~~
