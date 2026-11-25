# An endpoint reduction through $10^{20}$ and a critical-exponent bound

Let $E=(1,1,N-2)$. For each row write $a\le b\le c$ and $s=a+b$.

1. **Every counterexample at $N\le10^{20}$ must contain $E$.** In particular the full conjecture holds throughout this range whenever $N-1$ has at most two distinct prime divisors, by the [endpoint factor assignment](endpoint_reduction.md). The reduction uses published prime-gap computations.
2. **Any counterexample omitting $E$ has $N>10^{20}$, and each of its rows satisfies**
   $$s>10^8(N/10^{20})^{199/500}.$$
   In particular, each pair of smaller parts sums to more than one hundred million.
3. **For every sufficiently large $N$, each row other than $E$ in any counterexample satisfies**
   $$s>2.73\,N^{2/5}.$$
   More generally $2.73$ can be replaced by any constant strictly smaller than $C_*$ in (3) below. The asymptotic result is independent of the external prime-gap tables.

Consequently, for sufficiently large $N$, a prime power $Q=p^j\le N$ with $p\mid N(N-1)$ and $N-Q\le2.73N^{2/5}$ proves the full conjecture at $N$. A $p$-free row assigns the leading $Q$ to a single part since $N<2Q$. Its other two parts sum to at most $N-Q$, and $p$ divides the coefficient of $E$. Taking $p=2$ gives intervals above every sufficiently large power of two. This reaches the exponent excluded from the [earlier open range of exponents](power_intervals.md).

## Loss at small powers

Let $W(t)=t(t+1)/2$ and $B_t(x)=\prod_{j=0}^{t-1}(x-j)^{t-j}$. Give each row cutoffs $t_i$ summing to $H$ and fitting its parts, and write $P$ for the product of all nine $B_{t_i}(x_i)$. For $w=p^e$ set
$$F_w(d)=\sum_{j\ge0}(d-jw)_+,\qquad G_w(d)=F_{3w}(d).$$
Define $D(H)$ by assigning to each layer $w\le H$ the exponent
$$\max\bigl(0,\max_{0\le h<w}\{F_w(H-h)-G_w(H-h)-2G_w(H-m_w(h))\}\bigr),$$
where $m_w(h)$ is the greatest integer at most $3w-3$ congruent to $h$ modulo $w$. Multiply the corresponding powers of the prime bases.

For a counterexample,
$$ (N-H+1)^{W(H)}\le B_H(N)\le D(H)P. \tag{1}$$
Indeed a layer at $h=N\bmod w$ contributes $F_w(H-h)$ to the valuation on the left. A row avoiding $p$ contributes at least $G_w(H-h)$ to the valuation of $P$, and each other row contributes at least $G_w(H-m_w(h))$. At $w>H$ the avoiding row covers the entire layer. These lower bounds hold at every layer of each prime. [Layer proof](full_range_1b.md)

## Four inequalities and a prime-gap record

Suppose $E$ is absent. The [one-billion theorem](full_range_1b.md) gives $N>X=10^9$. By Bergman’s Proposition 10 every row has $s\ge11$. At most one row contains a part at most 2, since the full powers of $N$ assigned to two such rows would together divide 2. For a height divisible by 3 a row with minimum at least $H/3$ costs at most $(N/3)^{3W(H/3)}$. Balanced cutoffs on the two larger parts after using a fixed minimum give the analogous $(N/2)$ bound. [Bergman, Proposition 10](https://arxiv.org/pdf/0806.0607)

Assume a selected row has $s\le S=2000$. The table gives an upper bound for $P$ in (1).

| Selected minimum | $H$ | Selected cutoffs | Upper bound for $P$ |
|---|---:|---|---|
| $a=1$ | 12 | $(1,9,2)$ | $S^{45} N^{63}/3^{60}$ |
| $a=2$ | 6 | $(2,3,1)$ | $4 S^6 N^{19}/3^{18}$ |
| $a=3$ | 9 | $(3,5,1)$ | $108 S^{15} N^{39}/(3^{18}2^{20})$ |
| $a\ge4$ | 9 | $(4,4,1)$ | $(S/2)^{20} N^{39}/(3^{18}2^{20})$ |

At height 12 both other minima are at least 3. For a minimum equal to 3 use $(3,5,4)$ at cost $108(N/2)^{25}\le(N/3)^{30}$. At height 9 at least one other minimum is at least 3. The row with unrestricted minimum costs at most $(N/2)^{20}$, by taking $(1,4,4)$ or $(2,4,3)$ when the minimum is 1 or 2. For the last selected cutoffs use $ab\le(S/2)^2$. All cutoffs fit by $s\ge11$ and the ordering of the parts.

The constants $D(6),D(9),D(12)$ are $540$, $26{,}460{,}000$, and $107{,}820{,}266{,}400{,}000$. In each line $D(H)P<(X-H+1)^{W(H)}$ when the bound is evaluated at $N=X$. Its degree in $N$ is at most $W(H)$, so (1) fails at every $N\ge X$. Thus each such row has $s>2000$.

A row avoiding the largest prime $r\le N$ must have $s\le N-r$ by Bertrand’s postulate and Kummer’s criterion. The Prime Gap Search project reports exhaustive coverage through $10^{20}$. The largest consecutive-prime difference there is 1724. Thus $N-r<2000$ throughout that range, ruling out a selection without $E$. This use of the gap table is external to our one-billion enumeration. [Coverage records](https://primegap-list-project.github.io/fully-analyzed/) · [Maximal-gap table](https://pzktupel.de/RecordGaps/GAP01.php)

## The explicit bound above $10^{20}$

Here take $X=10^{20}$ and $S=10^8$. No part $m$ of a counterexample without $E$ can lie between 4 and 1000. At height 6 assign cutoff $\min(m,6)$ to it. If $m=4$ or $5$, put the remaining cutoffs on distinct other parts. The other two rows together cost at most $N^{18}/6^9$. Evaluating (1) at $X$ excludes the indicated interval, and the degrees preserve all the inequalities for larger $N$. More generally, any part $m\ge6$ satisfies
$$m\ge(6^9/540)^{1/21}(1-5/N)N^{1/7}.\tag{2}$$

Suppose $s\le S(N/X)^{199/500}$ in the selected row. At least one of the other rows has all its parts above 1000. Otherwise both have a part at most 3, and the product of their assigned powers of $N$ is at most 9. The selected parts would all be multiples of at least $N/9>s$. Call the row with three large parts $C$, and the other row $B$. In the table, $a$ is the selected minimum, $m$ is the minimum of $B$, and a star denotes a minimum above 1000. Distinct rows cannot both have minimum at most 2. All second parts exceed 1000 because their small-part sums exceed 2000.

| $a$ | $m$ | $H$ | Cutoff $k$ on $N-s$ |
|---|---|---:|---:|
| 1 | 3 | 6 | 0 |
| 1 | $*$ | 21 | 5 |
| 2 | 3 | 6 | 0 |
| 2 | $*$ | 9 | 2 |
| 3 | any | 6 | 0 |
| $*$ | 1 | 717 | 118 |
| $*$ | 2 or 3 | 6 | 0 |
| $*$ | $*$ | 9 | 1 |

Give $C$ equal cutoffs. If $m$ is small, give it cutoff $m$ and split the rest as equally as possible, with the larger cutoff on the smaller part. Otherwise use equal cutoffs. In the selected row with small $a$ take $(a,H-k-a,k)$. With large $a$ split $H-k$ equally up to rounding on the two smaller parts, putting the larger cutoff first. The corresponding bounds use $B_a(a)s^{W(H-k-a)}$ or $(s/2)^{W(u)+W(v)}$ on these two parts.

For each case, the checker verifies $D(H)P<(X-H+1)^{W(H)}$ at $N=X$, $s=S$, as an inequality on the indicated upper bounds. For a cost proportional to $N^d s^f$ it also verifies $500d+199f\le500W(H)$. These two tests prove the contradiction for every $N\ge X$. The large powers are compared by integer bit lengths giving upper and lower bounds on $1024\log_2 n$.

## The exponent $2/5$ with a constant

The key estimate when the height varies is
$$\limsup_{H\to\infty}\frac{\log D(H)}{H^2}\le d_0,
\qquad d_0=\frac12\sum_{j\ge0}\left(\frac1{3j+1}-\frac1{3j+2}\right)
=\frac{\pi}{6\sqrt3}.$$
To see this, use $m_w(h)\le h+2w$ and $F_w(x)=G_w(x)+G_w(x-w)+G_w(x-2w)$. Each layer loss is at most $G_w(H-w)-G_w(H-2w)$ by convexity of $G_w$. For fixed $k$, the prime number theorem in the form $\psi(y)=\sum_{p^e\le y}\log p\sim y$ gives
$$\sum_{w=p^e}(H-kw)_+\log p=k\int_0^{H/k}\psi(t)\,dt\sim H^2/(2k).$$
Use this at $k=3j+1$ and $3j+2$ for each fixed $j$. In the tail $j\ge J\ge1$, only $w\le H/(3J+1)$ occur. Each difference of positive parts is at most $w$, with at most $H/(3w)+1$ nonzero terms per $w$. Thus the tail is $O(H\psi(H/(3J+1)))=O(H^2/J)$. This justifies passage to the series. Integrating $1/(1+x+x^2)$ from 0 to 1 evaluates the sum before the factor $1/2$.

Suppose $s\le c_0N^{2/5}$ along an unbounded sequence of counterexamples omitting $E$, where $c_0$ is fixed. By (2), each part is at most 3 or grows at least as a positive constant times $N^{1/7}$. Classify $a,m$ as before. One of the ordinary rows has only large parts since $s<N/9$ eventually. Except when $a$ is large and $m=1$, fixed cutoffs already give total degree strictly below $W(H)$ at exponent $2/5$. Use the preceding table with three changes: for $(a,m)=(1,*),(2,*),(*,*)$, respectively, take $(H,k)=(24,5),(12,2),(12,0)$. The checker verifies the eleven strict degree comparisons.

In the remaining case take $H=12t$ with $t=\lceil(\log N)^2\rceil$. Put $(5t,5t,2t)$ on the selected row, $(1,6t,6t-1)$ on $B$, and $(4t,4t,4t)$ on $C$. These fit by (2) and the lower bound on second parts. Inequality (1) now gives
$$
(N-H+1)^{W(H)}
\le D(H)(s/2)^{25t^2+5t}N^{2t^2+t}
(N/2)^{36t^2}(N/3)^{24t^2+6t}.
$$
After substituting $s\le c_0 N^{2/5}$, the right side has excess degree $3t$ in $N$ over $W(H)=72t^2+6t$. Its contribution $3t\log N/H^2$ tends to zero. Replacing $N-H+1$ by $N$ on the logarithmic scale also costs $o(H^2)$. Taking logarithms and dividing by $H^2$ yields
$$0\le d_0+\frac{25}{144}\log c_0-\frac{61}{144}\log 2-\frac16\log3.$$
This excludes every
$$c_0<C_* =\exp\!\left(\frac{61\log2+24\log3-144d_0}{25}\right)>2.73.\tag{3}$$
When $E$ is present, each part of either other row exceeds $27^{1/5}(N-1)^{2/5}$. This gives the stated $2.73$ bound and also covers every $c_0<C_*$ since $C_*\approx2.73109<2\cdot27^{1/5}$. [Endpoint part bound](endpoint_factors.md)

Run `python3 check_small_part_bounds.py` to reproduce the finite inequalities, the degree checks, and a rational enclosure proving the strict margin in (3).
