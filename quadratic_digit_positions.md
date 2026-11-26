# Six powers at quadratic exponent positions

**Theorem.** For every prime $p$, every $a\ge1$, and all integers $t,c\ge0$, Wasserman's three-row conjecture holds at
$$ N=1+\sum_{k=0}^5 p^{a+kt+k^2c}. $$
The same common-prime conclusion holds for every choice of the three positive triples at each such $N$. Both coefficients of the quadratic and the outer shift are unrestricted. For example the theorem covers the six positions $a,a+t+c,a+2t+4c,a+3t+9c,a+4t+16c,a+5t+25c$ even when the position differences have gcd one.

The proof treats the two growth parameters together. At a root of the binomial for two selected positions, the magnitudes of the other four terms depend on the quadratic $k(k-i-j)$. Equal values of that quadratic pair indices by reflection. This leaves a small number of extra phase relations at 2 and 3, for which we give an exact finite verification.

## Reduction to a binomial and two short sums

If $t=0$ or $c=0$, the six positions lie in a progression of normalized span at most 25, or all coincide. These cases are covered by the [short-span](seven_span_and_gaps.md) and [three-position](seven_digit_units.md) theorems. Henceforth $t,c\ge1$. Put $h_k=kt+k^2c$ and $H=h_5=5t+25c$. All six positions are distinct. At $p\ge29$ the [effective-base theorem](seven_base27.md) already applies.

In the [binary allocation with 1 in the triple](seven_binary_alignment.md) the six positions split into gaps $n,n,n\pm1$. Each gap is $(j-i)(t+(i+j)c)$ for $i<j$. Equal gaps between two different pairs must have different coefficients of $t$: at the same width, equal sums of endpoints would give identical pairs. Their equality gives a nontrivial equation for $t/c$ with numerator of absolute value at most 25 and denominator at most 5. The two gaps differing by one give $\gcd(t,c)=1$. Thus $t\le25$, $c\le5$, and direct enumeration of the labelled pairs in this rectangle gives no such partition.

Every other possible counterexample has an unshifted pair on indices $i<j$, an unshifted trinomial on a three-element set $I$, and one remaining singleton. Write $R=\{0,\ldots,5\}\setminus\{i,j\}$, so $I\subset R$. Let $x,y,z$ be the sums of $p^{h_k}$ on the pair, trio and singleton. The [seven-unit allocation](seven_span_and_gaps.md#which-allocations-need-checking) requires $U\mid x,y$, $A\mid x,y+z$, $C\mid y,x+z$, all prime to $p$, with $U>1$ and $27(x+y+z)^2\le U^3A^2C^2$. The [two-trinomial bound](two_trinomial_cuts.md) gives $C\le4p^{2H/3}$. In particular
$$ 3\log_p U+2\log_p A>2H/3. \tag{1}$$

## A degree bound on a rational direction

Put $n=j-i$, $B=i+j$, $L=t+Bc$, and $b_k=k(k-B)$. Choose nonnegative coprime integers $q,r$ and put
$$ v=q+Br>0,\quad l=L/v,\quad \Delta=cq-tr,\quad \delta=\Delta/v=c-lr,\quad s_k=kq+k^2r=kv+b_k r. $$
Thus $h_k=l s_k+\delta b_k$. At $S=p^l$, the pair after stripping its lowest power is $F(S)=1+S^{nv}$. For $K=I$ or $R$ use $G_K(S)=\sum_{k\in K} p^{\delta b_k-m_K} S^{s_k}$, where $m_K=\min_{k\in K}\delta b_k$. These normalizations change the evaluated sums by powers of $p$ only. Let $d_K$ be the degree of the monic gcd of $F,G_K$ over $\mathbb Q(p^{1/v})$, and $w_K=\max_{k\in K}b_k-\min_{k\in K}b_k$. Then
$$ 3\log_p U+2\log_p A\le (3d_I+2d_R)L/v+9nv+E|\Delta|,\qquad E=n(3w_I+2w_R). \tag{2}$$

Here is a proof of (2) that also handles an irrational $p^l$. In the integer ring of this field, multiplication by $G_K$ modulo the monic polynomial $F$ has rank $nv-d_K$. Its rational image contains the monic gcd $h$ by Bézout. The coefficients of $h$ are algebraic integers since its roots are roots of unity. Choose independent matrix columns spanning the image and a nonsingular square minor $D$ of them. By Cramer's rule $Dh$ is in their span over the integer ring. For rank zero take $D=1$, $h=F$. At $S=p^l$ the assigned integer ($U$ or $A$) divides $Dh(p^l)$ after inverting $p$, and hence in the integer ring since it is coprime to $p$. The value is nonzero because $p^l>1$ and $h\mid F$.

For every complex embedding, each matrix column has Euclidean norm at most $|K|p^{|\delta|w_K}$, by signed reduction of its terms modulo $F$. Also $|h(p^l)|$ and every conjugate absolute value are at most $(1+p^l)^{d_K}$. Hadamard's bound and the field norm therefore give
$$ \log_p(U\ \text{or}\ A)\le d_K l+nv\log_p |K|+nv|\delta|w_K. $$
Here $d_K\log_p 2$ was included in $nv\log_p|K|$ using the rank $nv-d_K$. Summing gives (2), since $3\log_2 3+2\log_2 4<9$.

## Computing the common degrees

Suppose $\delta\ne0$. A common root of $F,G_K$ is a root of unity $\xi$ with $\xi^{nv}=-1$. Its terms have heights $\delta b_k$ and phases $\xi^{s_k}$. Two distinct indices have the same height exactly when their sum is $B$. No height occurs more than twice. If every height group vanishes separately, the trio gives no roots. For the four terms, all six indices (including $i,j$) would pair by reflection around $B/2$. This requires $B=5$. The three pair widths are $1,3,5$ and the slope differences on these pairs are their widths times $v$. The opposite-phase conditions together with $F=0$ are exactly $\xi^v=-1$. Thus the contributions from cancellation within heights are $d_I=0$ and $d_R=v\mathbf1_{B=5}$.

The [separation theorem](binomial_sparse_gcd.md) accounts for every root in this way at $p\ge5$. Any other root at 2 or 3 gives a minimal vanishing sum with unequal heights. For the trio, the sorted doubled heights after subtracting the minimum are $(0,0,1)$ or $(0,0,2)$ at 2, and $(0,0,1)$ at 3. The two lowest terms have relative phase with angle $\pm\pi/2$, zero, or $\pm\pi/3$, respectively, since their sum has magnitude $\sqrt2$, 2, or $\sqrt3$. Including the phase of the opposite top term, the relative phases all lie in the groups of orders $M_2=24$, $M_3=12$. For the four terms the [exact phase classification](quadrinomial_small_bases.md) gives the permitted height patterns and the same group orders.

For each choice of prime, pair and subset $K$, try both signs of $\delta$. The height range fixes its magnitude to the maximum listed doubled height divided by $2w_K$. Retain each match of the complete sorted list. Enumerate phase tuples in the indicated root group, anchoring the phase of the first index at one, and discard opposite pairs. Equality of each sum is tested modulo the cyclotomic polynomial, using $\sqrt p=\zeta_8+\zeta_8^{-1}$ or $\zeta_{12}+\zeta_{12}^{-1}$ as appropriate.

For a retained tuple, let $k_0$ be that first index. Write $a=(nv,(s_k-s_{k_0})_{k\ne k_0})$ and $b=(M_p/2,(b'_k)_{k\ne k_0})$, where $b'_k$ are the phase exponents modulo $M_p$. We need $a_h\theta=b_h/M_p\bmod1$ for every component. With $g=\gcd(a)$ there are exactly $g$ solutions if
$$ (a_h/g)b_j=(a_j/g)b_h\pmod{M_p}\quad\text{for all }h,j,$$
and none otherwise. Indeed an integral Bézout combination of the primitive vector $a/g$ solves for $g\theta$ modulo one. Distinct tuples give disjoint roots, and $F$ is squarefree. Adding these counts to the within-height contribution gives the exact degrees in (2).

## All sufficiently large spans

Choose a Dirichlet approximation $1\le v\le100$ with $|vc/L-r|\le1/100$. It has $r\ge0$ and $q=v-Br\ge0$, since $q\ge-B/100>-1$. Divide $q,r,v$ by the common content, which preserves the slope and the inequality. Then $|\Delta|\le L/100$. If $\Delta=0$, $(t,c)$ is an integer multiple of $(q,r)$, and $H=\gcd(t,c)(5q+25r)\le2500\gcd(t,c)$. At $H\ge68333$ the effective base $p^{\gcd(t,c)}$ exceeds 27.

When there are no additional phases, (2) has first term $2L\mathbf1_{B=5}=2H\mathbf1_{B=5}/5$. For the 15 pairs the verifier checks the rational inequalities
$$ 900n < 68333\left(2/3-(2/5)\mathbf1_{B=5}-(5nw_R/100)\max(1/5,B/25)\right). $$
Since $L/H\le\max(1/5,B/25)$ these bounds contradict (1) for $H\ge68333$.

For the possible extra phases, enumerate primitive $(q,r)$ with $v\le100$ and integral $v\delta$. Write $D=3d_I+2d_R$, $h=5q+25r$, $f=D/h$. The identity $H=lh+(25-5B)\delta$ makes the right side of (2) equal to
$$ fH+9nv+E|\Delta|-f(25-5B)\delta. $$
Among the 10,797 inputs tested, whenever $f<2/3$ the quotient of the constant term by $2/3-f$ is at most $136665/2$. The only cases with $f\ge2/3$ are in this table. Each has $r=1$, so $t=qc-\Delta$.

| $p$ | $(i,j)$ | $I$ | $q$ | $\Delta$ | $d_I,d_R$ |
|---|---|---|---:|---:|---|
| 2 | (1,3) | (0,4,5) | 1 | 1 | 10,0 |
| 3 | (1,4) | (0,3,5) | 7 | −1 | 12,12 |
| 3 | (1,4) | (2,3,5) | 7 | 1 | 12,12 |

For each row, form the two trinomials for $y,x+z$ in $T=p^c$ with terms $p^{-k\Delta}T^{kq+k^2}$ and clear monomial powers of $T$ and $p$ separately. The verifier finds polynomial gcd one in each case. Both degrees are at most 60, and each polynomial has coefficient 1-norm at most $3p^5$, since $|\Delta|=1$ and $0\le k\le5$. Their nonzero resultant gives $\log_p C\le120\log_p(3p^5)<840$. Also $f=1$ in all three cases, and (2) gives $3\log_p U+2\log_p A\le H+1000$ there using $v\le12$, $E\le180$ and $|\delta|\le1$. The total is below $2H$ at $H\ge68333$. This finishes the large-span proof.

## Exact finite verification

For $5t+25c\le68333$ the verifier first applies (2) to intervals of $t$ at each fixed $c,i,j$. Use $E_0=5nw_R$ in place of $E$ to cover every singleton choice. For each direction with $v\le100$ it considers the interval specified by
$$ 3(9nv+E_0|cq-tr|)+6L\mathbf1_{B=5}\le2H. $$
It skips points with $cq=tr$. At the 72 direction and nonzero determinant entries where the degree computation can increase $D$ above $2v\mathbf1_{B=5}$, it uses the greatest computed $D$ over the two primes and all singletons. A point is retained in the interval only when $3(DL/v+9nv+E_0|cq-tr|)\le2H$. The C++ verifier evaluates that integer inequality individually at every point it marks. Intervals only accelerate the search for a witness. Every marked point contradicts (1).

Ignoring $\gcd(t,c)\ge5$, where the effective-base theorem applies, the interval calculation leaves 411,088 tuples $(i,j,t,c)$. The common step of the positions is $g=\gcd(t+c,2c)$. Apply the short-span theorem if $H/g\le100$ and the effective-base theorem if $p^g\ge27$. This leaves 2,216,569 base comparisons. At each remaining prime the integer verifier checks the assigned divisors directly. It reduces the tri- and quadrinomial modulo $1+p^{nL}$, rotates past a largest exponent gap by a unit power of $p$, and computes their gcds with that binomial. A one-term signed remainder also gives its absolute coefficient as an upper bound. A zero quadrinomial remainder forces $U=1$. Otherwise the verifier compares $16^3 U^9 A^6$ with $27^3p^{2H}$ using upper bounds on $U,A$ or their integer gcds. When necessary it then uses the full numerical $\gcd(y,x+z)$ and $m=x+y+z$ to test $U^3 A^2 C^2<27m^2$. The full product comparison is used 2,279 times. This completes the proof.

Run the following with Python, python-flint, SymPy and a C++17 compiler. It recomputes the phase table entries, checks all 10,797 pairs of degrees by independent polynomial gcds, compiles the interval enumeration, and streams the remaining tuples into the integer verifier. Separate direct tests recertify 1,200 sampled exclusions and check 983 full products using the original, unrotated inputs on sampled residual allocations.

~~~sh
python3 check_quadratic_digit_positions.py
~~~

## A minimum four-prime certificate at a 52-digit weight

Take $p=2$, $a=6$, $t=13$ and $c=4$. Here $m=1+2^{17}+2^{42}+2^{75}+2^{116}+2^{165}$ and $N=1+64m$. Four primes certify every three-row choice at this one weight:
$$2,\qquad Q=N/96537121,\qquad s=m/[3(2^{33}+1)],\qquad r=N-2.$$
The [integer certificate](check_quadratic_prime_list.py) proves primality by Pocklington's criterion with witness 2 throughout and trial division at the leaves.

Any row avoiding 2 assigns the top bit $2^{171}$ to one part. Its two smallest parts then sum to at most $N-2^{171}<\min(2Q,s)$. Avoiding $Q$ requires all three parts to be multiples of $Q$, and avoiding $s$ requires at least two to be multiples of $s$. So this row carries at both primes. A row avoiding both $Q$ and $s$ would have two positive multiples of $Qs$ and a third positive multiple of $Q$, whereas $N<Q(2s+1)$. Finally, the only row avoiding $r$ is $(1,1,r)$ up to permutation. It carries at 2, $Q$ and $s$. Thus each row avoids at most one of the four primes.

Four is minimal. Every prime $w$ admits a positive row avoiding it when the base-$w$ digit sum of $N$ is at least three, by distributing its digit units among three nonempty parts. Here $N$ has seven binary ones. An odd-base digit sum of two is impossible by parity, and $N$ is not a prime power since $Q\mid N$ and $1<N/Q<Q$. Given any list of at most three primes, one can therefore choose one row per prime to avoid it.

~~~sh
python3 check_quadratic_prime_list.py
~~~
