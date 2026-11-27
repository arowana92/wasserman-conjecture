# Three binary blocks with two independently varying gaps

Write
$$
 N=D2^a+m2^b+C,\qquad a\ge b\ge0,
$$
with positive integer coefficients. The blocks may overlap.

**Uniformity in the two outer coefficients.** For all but $O(X^{39/40})$ integers $m\le X$, Wasserman's three-row conjecture holds at **every** such weight for
$$
             1\le D\le X^{1/32},\qquad 2\le C\le X^{1/32}.
$$
Thus the same middle coefficient gives entire intervals of consecutive good weights, at every pair of exponents and every indicated prefix.

**Three unrestricted coefficients.** Among $1\le D,m,C\le X$, there are at most **$O(X^{2.98})$ ordered triples** for which even one exponent pair gives a counterexample. This counts the union over both unbounded exponents.

The key calculation for very large exponents gives a stronger bound at each fixed $D,C$ with $C\ge2$. Put $n=\lceil\log_2(2X+2)\rceil$ and take $X\ge\max(D,C)$. At most $O_{D,C}(X^\alpha\log(2X))$ middle coefficients can fail at any $a>512n$, where $\alpha=0.551463\ldots$ is the [binary affine-triple exponent](binary_affine_counts.md).

## The binary row

In a counterexample use the [endpoint reduction](all_weights_endpoint_reduction.md) and choose a nonendpoint row avoiding 2. Assign it the full prime powers with products $U\mid N$, $S\mid N-1$, assigning 2 to it. Every one of its parts $y$ is at least $U\ge3$ and satisfies $U\mid y$, $S\mid y(y-1)$. Split $S=S_1S_2S_3$ according to the coordinate having residue one. The [complementary-row size bound](projected_rows.md) gives
$$
                     27(N-1)^2\le U^3 S^2. \tag{1}
$$
For $Q=2^a,R=2^b$ the binary row can be written
$$
                    A_i=d_i Q+m_i R+c_i.
$$
The three columns separately split the bits of $D,m,C$ without a carry. To see this even when blocks overlap, merge pairs of equal powers of two until the binary units of $N$ remain. Assign the original units to the parts receiving their merged units in the carry-free row.

## A large divisor for a matrix of three columns

Assume $C\ge2$, $a>512n$, and both $b,a-b>4n$. By (1) one of $U,S_1,S_2,S_3$ is at least $(27(N-1)^2)^{1/9}$. Pair that factor with the vector $e=0$ or the corresponding coordinate unit vector, and write $\delta=\sum e_i\in\{0,1\}$. Let $W$ be its odd part. The power of two discarded is at most $2^n$ since $N\equiv C\bmod 2^b$ and $2\le C<2^b$. Thus
$$
 W\ge 2^{2a/9-n},\qquad
 W\mid d_iQ+m_iR+c_i-e_i\quad(1\le i\le3). \tag{2}
$$

Let the columns of the matrix be $d_i,m_i,c_i-e_i$. If it has rank one, there is a primitive positive direction $z$ and the three original columns are
$$
                         g z,\quad h z,\quad l z+e. \tag{3}
$$
Indeed the first column is nonnegative and nonzero, and any zero coordinate of its primitive direction would give an actual part at most one. All column multipliers are integers by primitivity, and an unchanged coordinate shows $l\ge0$.

Otherwise the first two columns must be independent. If proportional, eliminate both variable terms between two rows using the first column. Then $W$ divides every minor of columns one and three, one of them nonzero and of absolute value at most $2X^2$. This contradicts (2). With the first two columns independent, the determinant of the three-column matrix is divisible by $W$ and bounded in absolute value by $6X^3$. It vanishes. Cramer's rule gives a rational solution $(x,y)=(I/\Delta,J/\Delta)$ of all three linear equations in place of $(Q,R)$, with
$$
 0<|\Delta|\le2X^2,\quad |I|,|J|\le2X^2,\quad
 W\mid\Delta Q-I,\quad W\mid\Delta R-J.
$$
Neither numerator vanishes since $W$ is odd and greater than $|\Delta|$.

Dirichlet's pigeonhole argument applied to $b/a$ gives $1\le u\le32$, $0\le v\le u$ and $|r|=|v a-u b|\le a/32$. Compare $I^v\Delta^{u-v}$ with $2^r J^u$, clearing a negative power of two if needed. Their integer difference after clearing is divisible by $W$ and has absolute value at most $2^{64n+1+a/32}<W$. Thus $x^v=2^r y^u$ and $|r|\le 2nu$. The cases $v=0$ and $v=u$ would force $b\le2n$ or $a-b\le2n$, respectively. Reduce $u,v$ by their gcd in the relation on **absolute values**. The odd-prime valuations then give
$$
 x=\pm2^j t^u,\qquad y=\pm2^k t^v,\qquad 0<v<u\le32, \tag{4}
$$
where $t>0$ is rational with odd numerator and denominator and $|j|,|k|\le 2n$. In this formula $u,v$ are the reduced indices.

Summing the equations gives $Dx+my+C-\delta=0$. Substitute (4) and clear powers of two. The rational-root theorem shows that the odd numerator of $t$ divides $C-\delta$ and its odd denominator divides $D$. For fixed $D,C$ there are therefore only
$$
 O\big(n^2\tau(D)(\tau(C)+\tau(C-1))\big)
$$
possible $m$ in the rank-two case. Here $\tau$ is the divisor-counting function. The middle coefficient is determined once the signs, indices, powers and $t$ are chosen. The assumption $C\ge2$ ensures a nonzero constant coefficient throughout.

For an example where a rank-two divisor really grows, take the rows of the coefficient matrix to be $(1,1,18)$, $2^7(1,2,9)$ and $2^{14}(2,1,45)$, with $e=0$. At exponents $(3k,2k)$ the three parts have common divisor $2^k+3$, since the common point is $(-27,9)$. The columns split their binary digits freely, and so do the actual parts for $k\ge24$. This example concerns the linear congruences in (2).

## When one gap is short

Still suppose $a>512n$. If $b\le4n$, merge the bottom two blocks to get $D2^a+(m2^b+C)$. If $a-b\le4n$, merge the top two blocks to get $(D2^{a-b}+m)2^b+C$. In each case the two resulting blocks are separated, and their product is at most $2^{6n}$. Its fifth power is less than $27(N-1)^2$. The [two-block determinant theorem](projected_rows.md#two-blocks-of-digits-in-any-prime-base) forces the upper parts to be $g z$ and the lower parts to be $l z+e$, with $z$ primitive positive and $e$ zero or one coordinate unit. Both triples add without a binary carry. The sum $s$ of $z$ is odd, since its coordinates have distinct 2-adic valuations and one is odd.

First let $D,C\le Y=X^{1/32}$ and $C\ge2$. In the first merger $s\mid D$; in the second $s\mid C-\delta$ and $C-\delta>0$. Thus $s\le Y$, giving $O(Y^3)$ directions. In the first merger $sl=m2^b+C-\delta$, so for fixed $b,C,s$ the possible $l$ run in one progression of step $2^b$ with at most $1+X/s$ terms. In the second, $sg=D2^{a-b}+m$, so $g$ is in an interval of at most $1+X/s$ terms for fixed $a-b,D,s$.

On either progression only $O(X^\alpha)$ inputs can make the varying triple add freely. In the [affine-triple count](binary_affine_counts.md), substitute the initial term plus the step times an index. The spread of the 2-adic valuations of the slopes, including their sum, is still at most $\log_2 s$, and at most $\log_2(1+X/s)+O(1)$ index bits are needed. Thus all short gaps together contribute $O(Y^4 X^\alpha\log(2X))$ middle coefficients. For the rank-one case (3), $s\mid D$ and $m=hs$ with $hz$ carry-free, giving $O(Y^3X^\alpha)$. The rank-two alternatives (4) contribute $Y^2X^{o(1)}$. With fixed $D,C$ these estimates also give the claimed $O_{D,C}(X^\alpha\log(2X))$.

For the unrestricted box discard $C=1$, costing at most $X^2$ triples. Use $Y_0=X^{1/10}$. For $s\le Y_0$, count just the end coefficient given by $gz$ or $lz+e$ in the merged shape, or by $gz$ in (3). It has $O(X^\alpha)$ possibilities at a fixed direction. This costs $O(Y_0^3 X^{2+\alpha}\log(2X))$. For $s>Y_0$ one has either $s\mid D,\ s\mid m2^b+C-\delta$, or $s\mid C-\delta,\ s\mid D2^{a-b}+m$. At each $s\le X$ and fixed short gap there are $O((X/s)X(1+X/s))$ triples, by fixing the unmapped coefficient and using residue classes for the other two. The same estimate holds directly from $s\mid D,m$ in (3). Summation gives $O((X^3/Y_0+X^2\log(2X))\log(2X))$. The rational-root count contributes $X^{2+o(1)}$. This proves the unrestricted estimate for large heights.

## Bounded heights with small outer coefficients

Now $a\le512n$, so there are $O(n^2)$ exponent pairs. All integers being factored have polynomial size in $X$. For $K=US\le X^{29/32}$, the [two-congruence reconstruction applied to the other row](two_predecessor_congruences.md#recovering-the-row-from-a-product) gives $X^{29/32+o(1)}$ possible weights. Apply that bound on dyadic intervals, counting weights below the cutoff directly. At fixed exponents, $D,C$ each weight gives at most one $m$. Allowing $D,C\le Y$ thus costs $X^{31/32+o(1)}$ altogether.

Write $H(t)=-t\log_2 t-(1-t)\log_2(1-t)$. At most $O(X^{H(3/5)})$ middle coefficients have more than $3n/5$ ones. In the other cases one $m_i$ has at most $n/5$ ones. The binomial tail bound gives $O(X^{H(1/5)})$ such values. Its part determines $(U,S)$ up to $X^{o(1)}$ divisor choices. For fixed $D,C$ that part has at most $(D+1)(C+1)$ choices of its two outer contributions. Each label fixes $N\bmod K$ by CRT. Since $\gcd(2^b,K)\le C$ (the shared power divides either $C$ or $C-1$), there are at most $1+X C/K$ possibilities for $m$. The total in the large-$K$ case is
$$
 X^{H(3/5)+o(1)}+X^{H(1/5)+o(1)}Y^4(1+XY/X^{29/32}).
$$
Both exponents are strictly below $39/40$. This proves the uniform middle-coefficient assertion.

## Bounded heights in the unrestricted box

For fixed $a,b$ carry the original coefficients into three nonoverlapping fields: reduce $C$ modulo $2^b$, add the quotient to $m$, then carry its quotient modulo $2^{a-b}$ into $D$. The three field widths as counts of possible values are bounded by the powers of two
$$
 L_0=2^{\min(b,n)},\qquad L_1=2^{\min(a-b,n+1)},\qquad L_2=2^{n+2}.
$$
Put $V=L_0L_1L_2\gg X$. A weight has $O(X^3/V)$ original representations, by fixing $C$ in a class modulo $2^b$ and $m$ subsequently in a class modulo $2^{a-b}$. The possible bit positions form three disjoint fields of $\log_2 V$ bits in total. A field of width one is empty.

Discard original tuples with $2^{v_2(N(N-1))}>X^{1/40}$. There are $O(X^{3-1/40}+X^2)$, by fixing $D,m$ and counting $C$ in two residue classes. Write $K_o$ for the odd part of $K$. When $K_o\le V X^{-1/20}$, the reconstruction bound used above gives at most $VX^{-1/40+o(1)}$ weights. For larger $K_o$ use the lattice of integer vectors $z$ satisfying
$$
                  z_0+2^b z_1+2^a z_2\equiv0\bmod K_o,
 \qquad \|z\|=\max_i |z_i|/L_i.
$$
Each label restricts the encoded weights to one affine class. Write $\lambda$ for the shortest nonzero vector length.

When $\lambda>X^{-3/40}$ a packing into small boxes gives $O(X^{9/40})$ points per translate. Outside $O(V^{H(3/5)})$ weights, one part has at most $(\log_2 V)/5$ binary ones in these fields, hence at most $O(V^{H(1/5)})$ values. Each has $X^{o(1)}$ labels. Since $V\gg X$, the count is $O(VX^{-1/40+o(1)})$.

For $\lambda\le X^{-3/40}$ take dyadic $r\le\lambda<2r$ and let $I$ be the coordinates with $2rL_i\ge1$. Its size $k$ is at least one. A shortest vector is supported on $I$, with $O(r^k V_I)$ choices, $V_I=\prod_{i\in I} L_i$. Its dot product with $(1,2^b,2^a)$ is nonzero: reducing successively modulo $2^b$ and $2^{a-b}$ forces the first two coordinates to vanish if the dot product vanishes, since $2r<1$ for large $X$. Thus there are at most $X^{o(1)}$ choices for its divisor $K_o$ and the label.

Fix the coordinates outside $I$. Each translate in the remaining box has
$$
                         O(r^{-(k-1)}+V_I/K_o)
$$
points. In affine dimension below $k$, pack points at separation at least $r$ inside the bounded scaled box in that affine subspace. In full dimension fix one point, and choose $k$ differences to it having maximal absolute determinant. Every other difference has coordinates of absolute value at most one in that basis by Cramer's rule. The difference lattice in the slice has index $K_o$ since each dot-product coefficient is a unit modulo the odd $K_o$. The chosen differences thus span a sublattice of relative index at most $k!V_I/K_o$. Each coset supplies at most $3^k$ points in the indicated basis box. This proves the bound in all dimensions, including $k=1$.

Multiplying by the vector count and the $V/V_I$ outside choices gives
$$
               VX^{o(1)}(r+r^k V_I/K_o)\ \ll\ VX^{-1/40+o(1)}.
$$
Indeed $K_o>V X^{-1/20}$ and $r\le X^{-3/40}$. Sum the logarithmically many scales and exponent pairs and multiply by the representation bound. The resulting $X^{3-1/40+o(1)}$ is $O(X^{2.98})$. This completes the proof.

Run the standard-library verifier for the entropy margins, the lattice estimates and the rational relation, including a growing rank-two example.

~~~sh
python3 check_three_binary_blocks.py
~~~
