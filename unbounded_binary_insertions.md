# Inserting two binary gaps of unrestricted lengths

**Theorem.** Apart from $O(2^{0.98n})$ integers $w$ of binary length $n$, the following holds. Choose any two distinct cuts in the binary expansion of $w$ for which the lowest block has value at least two. Insert any numbers of zeros at the two cuts. Every resulting weight satisfies Wasserman's three-row conjecture. Either gap may be empty and the gap lengths have no upper bound. The statement also holds for one such cut or no cuts.

Thus one good input certifies a doubly infinite family for each eligible pair of cuts. For odd inputs the lowest block must retain at least two ones. For even inputs it is enough that the lowest block retain the least significant one.

Fix the cut positions. Write $D,m,C$ for the values of the high, middle and low blocks, and $D_0,M_0,C_0$ for two raised to their respective bit lengths. Then $1\le D<D_0$, $0\le m<M_0$, $2\le C<C_0$, and $V=D_0 M_0 C_0=2^n$. The output is
$$ N=D2^a+m2^b+C,\qquad b\ge\log_2 C_0,\quad a-b\ge\log_2 M_0. $$
We count the input triples allowing even one bad exponent pair. The proof treats unequal block lengths uniformly.

## A uniform bound on polynomial graphs

For any rational polynomial $P$ of degree between 2 and 32, the number of integer pairs $(t,P(t))$ in a square of side $L\ge1$ is $O(1+L^{2/3})$, with an absolute implied constant independent of the coefficients. Partition the horizontal range at real roots of $P'$ and $P''$, and count the endpoints separately. On each remaining monotone strictly convex or concave piece the consecutive lattice-point differences, ordered by increasing first coordinate, have distinct slopes. Their total $\ell^1$ length is at most $2L$ by monotonicity. At most $O(R^2)$ distinct integer differences have length at most $R$. At most $2L/R$ have greater length. Take $R=L^{1/3}$. This elementary convexity bound will allow the largest *outer* block to vary.

## Heights at most $512(n+2)$

All exponent pairs in this range give at most $O(n^2)$ fixed placements, each with disjoint free binary fields. Remove the $O(V^{39/40})$ input words divisible by a power of two greater than $V^{1/40}$. Every remaining even output has this same bound on its power of two, since $C\ge2$ and is even when the output is even. In fact the lowest one stays in the low block. For odd outputs use the unrestricted odd case of the [field-counting argument](binary_digit_insertions.md).

Specifically, in that argument take $r=3$, no fixed ones, $h\le514(n+2)$, $Y=V^{4/5}$ and $\epsilon=1/40$. All the required bounds still hold:
$$
H(3/5)<39/40,\qquad H(1/5)+9/40<39/40,\qquad
8/5-(2/3)(39/40)<39/40.
$$
The volume term uses $H(1/5)+1/5<39/40$ and the odd divisor-function losses are $V^{o(1)}$. This counts at most $V^{39/40+o(1)}$ possible input words for all bounded placements.

## Very large heights

Consider $a>512(n+2)$ and first suppose $m>0$. The equations below are consequences of the [three-block rank calculation](three_binary_blocks.md), applied with $n+2$ as an upper bound on its logarithmic coefficient parameter. Its estimates need only the upper bound $2^n$ on the coefficients.

When both $b,a-b>4(n+2)$ the binary carry-free row has digit columns of one of two types. In rank one these are $gz,hz,lz+e$, with $z$ a positive primitive triple and $e$ zero or a unit coordinate vector. In rank two there exist signs $\sigma,\tau$, relatively prime integers $0<v<u\le32$, integers $|j|,|k|\le2(n+2)$ and a positive reduced odd fraction $t=t_1/t_2$ such that
$$
 D\sigma2^j t^u + m\tau2^k t^v+C-\delta=0,\qquad
 t_1\mid C-\delta,\quad t_2\mid D,\qquad \delta\in\{0,1\}. \tag{1}
$$
All these alternatives hold uniformly over both exponents. If either gap is at most $4(n+2)$, merge its two blocks. The two resulting digit columns of the row are $gz,lz+e$. There are only $O(n)$ short-gap parameters, and no dependence of this shape on the length of the other gap. Again $z$ is positive and primitive. Put $s=\sum z_i$ and $\delta=\sum e_i$. The triple $gz$ adds freely in binary, hence its three nonzero coordinates have distinct 2-adic valuations. The primitive direction has exactly one odd coordinate and $s$ is odd.

**Rank two with arbitrary block widths.** Put $L=\max(D_0,M_0,C_0)\ge V^{1/3}$. If $L=M_0$, fix $D,C$. By (1) there are $V^{o(1)}$ choices of the fraction and the other displayed parameters. Each choice gives at most one $m$. If $L=D_0$, fix $m,C$ and the numerator $t_1\mid C-\delta$. The denominator $1\le t_2\le D<L$ varies, and
$$D=-\frac{(C-\delta)t_2^u +m\tau2^k t_1^v t_2^{u-v}}{\sigma2^j t_1^u}.$$
This is a rational polynomial in $t_2$ of degree $u\ge2$, since $C-\delta>0$. The graph bound gives $O(1+L^{2/3})$ values of $D$ for each choice. If $L=C_0$, fix $D,m$ and $t_2\mid D$ instead. Equation (1) makes $C$ a degree-$u$ rational polynomial in $1\le t_1\le C<L$, with nonzero leading coefficient. The same bound applies. Thus all rank-two cases together cost at most $V^{1+o(1)}L^{-1/3}\le V^{8/9+o(1)}$ input words.

**Rank one and short gaps.** Split at $s=V^{1/32}$. Below that cutoff there are $O(V^{3/32})$ choices of $z$. Fix the two coefficients outside a block of largest width $L$, and the short-gap parameter if present. Its varying merged or unmerged coefficient has shape $sg$ or $sl+\delta$ and is an affine function of the original block with slope 1 or a power of two. Since $s$ is odd, the varying multiplier $g$ or $l$ is restricted to a progression of at most $1+L/s$ terms. Only $O(L^\alpha)$ values can make the corresponding triple add freely, where $\alpha=0.551463\ldots$. Indeed after substituting the progression into $g z$ or $l z+e$, the 2-adic spread of the slopes including their sum is at most $\log_2 s$, and at most $\log_2(1+L/s)+O(1)$ index bits are needed. Apply the [affine-triple prefix bound](binary_affine_counts.md). Here $s\le V^{1/32}\le L$. Multiplying the counts gives at most
$$ V^{1+3/32+o(1)}L^{\alpha-1}\le V^{1+3/32-(1-\alpha)/3+o(1)}. $$

For large $s$, in a bottom merger $s\mid D$ and $s\mid m2^b+C-\delta$. For each fixed $m,C,b,\delta$ the second integer is positive and has $O(n)$ bits. It thus allows $V^{o(1)}$ possible divisors $s$. Each gives at most $D_0/s$ positive multiples for $D$. In a top merger $s\mid C-\delta$ and $s\mid D2^{a-b}+m$. Fix $D,m,a-b,\delta$ and count $C$ using $C-\delta>0$. In the three-column rank-one case use $s\mid D,m$ directly. The resulting total is $V^{31/32+o(1)}$. These divisor counts avoid any assumption that the two narrower blocks have comparable widths.

If $m=0$, the two gaps concatenate with the zero middle block into one. For a single cut the two-block determinant theorem at large height gives $D=g s$, $C=l s+\delta$. The same argument with $V'=D_0 C_0$ instead of $V$ counts the large heights. In the small-$s$ count $L\ge(V')^{1/2}$, and in the large-$s$ count fix $C$ and divide $C-\delta$. If $V'\le V^{39/40}$ one may simply count the $V'$ input choices. Otherwise the polynomial and divisor losses are still $(V')^{o(1)}$. This also handles inputs originally given with a single cut. For no cuts the bounded-height field estimate applies to one unmodified field.

All the large-height exponents are less than $39/40$. There are only $O(n^2)$ cut pairs, so the final bound $V^{39/40+o(1)}$ implies the theorem.

Run the integer tests for the two polynomial parametrizations, the convexity count and the exponent margins.

~~~sh
python3 check_unbounded_insertions.py
~~~
