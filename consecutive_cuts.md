# Consecutive cuts in the two nonendpoint rows

Suppose $2\le a<a+1<N$ and
$$\gcd\left(N(N-1),\binom Na,\binom N{a+1}\right)=1.$$
The following are necessary:
$$m=\frac{a(a+1)}N\in\mathbb Z,\quad m>10^6,\quad
q=\frac{m(m-2)}{N-1}\in\mathbb Z,\quad q\ge13. \tag{1}$$
Consequently $m\ge1+\sqrt{1+13(N-1)}$. The million bound holds at every weight, and the lower bound on $q$ holds for arbitrarily large $m$. As $N$ grows one also has the effective bound
$$a\gg N^{3/4}\frac{(\log N)^{1/36}}{(\log\log N)^{1/4}}. \tag{2}$$
Each conclusion applies to any two consecutive cuts in a possible trinomial counterexample. Either cut may be a single part or the sum of two parts.

## The divisor and the common gcd

At every full power of a prime of $N$ one cut must be zero. Hence $N\mid a(a+1)$. At every full prime power of $N-1$, either $a$ or $a+1$ must have residue 0 or 1. Thus $a$ has residue $-1,0$ or 1, and $m\equiv a(a+1)$ has residue 0 or 2. This proves $N-1\mid m(m-2)$. The [quotient bound](cross_cut_quotients.md) requires $m>1$. If $m=2$, the integer $k=m(N-a)(N-a-1)/(N(N-1))$ lies strictly between 0 and 2. Substitution of $k=1$ gives $N=4a-3$. Combining with $2N=a(a+1)$ yields $(a-1)(a-6)=0$. At $a=6,N=21$ both cuts carry in base two. Hence $q$ is positive.

For arbitrary consecutive cuts let $D,J$ be the gcds of the two binomials with $N,N-1$ respectively and let $v=a(a+1)/N=t/d$ in lowest terms. For $v\ne2$,
$$DJ\ge \sqrt{\frac{N-1}{|v(v-2)|}}. \tag{3}$$
Indeed $d\mid D$ since $a,a+1$ are coprime and $d$ is the common factor of the two reduced cut denominators. At each retained power in $(N-1)/J$, the same residue argument applies, so that power divides $t(t-2d)$. Use $(DJ)^2\ge d^2J$.

In the exceptional algebraic case $v=2$ there is the explicit divisor
$$\frac{a+2}{\gcd(a+2,12)}\mid J.$$
Here $N-1=(a-1)(a+2)/2$. At every layer $p^j\mid a+2$ with $p^j>3$ and $p^j\mid N-1$, the cut residues are $p^j-2,p^j-1$, both larger than 1. For odd $p$ all layers through $v_p(a+2)$ divide $N-1$; only the layer 3 may fail the size test. For $p=2$ one loses at most the top layer and the layer 2. Summing the remaining carries gives the divisor. In particular for any $0<\epsilon\le1/4$, two consecutive cuts with $a\le N^{3/4-\epsilon}$ have common endpoint gcd $\gg N^{2\epsilon}$, by this divisor and (3).

## An elliptic equation

The triangular-number identity and the definition of $q$ give
$$q(2a+1)^2=4m^3-8m^2+4qm+q. \tag{4}$$
The cubic on the right has discriminant $-16q(64q^2+107q-128)\ne0$ for every $q\ge1$. Write $q=dt^2$ with $d$ squarefree. The point $x=dm$, $y=d^2ta$ is integral on
$$y^2+d^2ty=x^3-2dx^2+qd^2x. \tag{5}$$

For an effective estimate put $c_4=d^2(64-48q)$, $c_6=d^3(512-792q)$. The integers $X=36x-24d$ and $Y=108(2y+d^2t)$ satisfy $Y^2=X^3-27c_4X-54c_6$. The polynomial discriminant has absolute value $4\cdot27^3\cdot1728\,d^6 q(64q^2+107q-128)$. Hajdu and Herendi's Theorem 1, applied to this nonsingular integral short equation, gives
$$\log m\ll d^6q^3(\log(2q))^9.$$
Indeed their $c_1$ is $O(|\Delta|^{1/2}\log^4(2|\Delta|))$, their $c_2$ has polynomial size in $q$, and their upper bound for the logarithmic size of a point is $O(c_1\log(c_1)(c_1+\log(c_2)))$. In particular $q\gg(\log m)^{1/9}/\log\log m$ since $d\le q$. Equation (2) follows by $m=1+\sqrt{1+q(N-1)}$ and $a(a+1)=mN$. For any fixed upper bound on $q$, the resulting bound on $m$ is effective. [Hajdu and Herendi, Theorem 1](https://math.unideb.hu/sites/default/files/inline-files/ellracuj.pdf)

There is an elementary improvement if $m=dv^2$ lies in the same square class as $q=dt^2$. Then $q>\sqrt{dm}$ and $m>d^{1/3}N^{2/3}$. To see this, (4) gives
$$[t(2a+1)]^2-[2v(m-1)]^2=4(q-1)v^2+t^2>0.$$
The gap is at least $4v(m-1)+1$. If $q\le dv$, then $t^2\le v$ and the right side is at most $4(dv-1)v^2+v$, smaller than the gap for $v>1$. For $v=1$ one has $q<m-2=d-2$ since $N-1>m$, which is incompatible with $q=dt^2$. Finally $N=1+m(m-2)/q<m^{3/2}/\sqrt d$.

## Fixed slopes with an offset of one

For every fixed integer $s\ge1$, the bound (2), with an implied constant depending on $s$, holds whenever the two endpoint-coprime cuts are $a$ and $b=sa+1$ and $N$ is sufficiently large. This uses a family of nonsingular cubics for every slope.

Here $m=ab/N$ is an integer since $\gcd(a,b)=1$. At a full predecessor power, if a cut is zero then $m$ is zero there. If $a$ has residue 1 then $m$ has residue $s+1$. If $b$ has residue 1 then $sa=b-1$ is zero, and hence $sm$ is zero. Thus
$$q=\frac{s m(m-s-1)}{N-1}\in\mathbb Z.$$
The [seventh-power bound](cross_cut_quotients.md) ensures $m>s+1$ for large $N$. Eliminating $N$ gives
$$q(2sa+1)^2=4s^2m^2(m-s-1)+q(4sm+1). \tag{6}$$

A repeated root $z$ of the cubic would be rational, since a repeated irrational root of a rational cubic would require its conjugates to be repeated as well. Neither zero nor $-1/(4s)$ can be a root. Solving (6) at a zero of its right side expresses $q=4s^2z^2(s+1-z)/(4sz+1)$. At a repeated zero the derivative of this rational function vanishes, giving
$$8s z^2-(4s(s+1)-3)z-2(s+1)=0.$$
Its discriminant is $H^2-16$ with $H=4s(s+1)+5\ge13$, strictly between consecutive squares. This excludes every repeated root.

For detail on uniformity in $q$, write the polynomial as $Az^3+Bz^2+Cz+D$ with $A,B$ fixed at slope $s$ and $C=4sq,D=q$. For $q=dt^2$ as above the integral substitution $Z=Ad m$, $V=Ad^2t(2sa+1)$ gives
$$V^2=Z^3+Bd Z^2+ACd^2 Z+A^2Dd^3.$$
Removing the quadratic term by $X=9Z+3Bd,\ Y=27V$ gives an integral short model. Its cubic discriminant is $O_s(d^6q^3)$ by the discriminant formula, and its coefficients have polynomial size in $q$. The same height bound therefore gives $\log m\ll_s q^9\log^9(2q)$. Since $m\gg_s\sqrt{qN}$ and $a(sa+1)=mN$, this proves (2) at the stated slope.

## Unbounded quotients at 27 values of the parameter

We used the published integral-point tables for elliptic curves over $\mathbb Q$ and compared the actual integral Weierstrass coefficients. Among $1\le q\le1000$, the following 27 values match a model in the tables after the integral transformations below:
$$1,2,3,4,5,6,7,8,9,10,11,12,14,16,18,20,24,25,27,32,36,40,56,64,72,128,384.$$
For these values all integral points yield just 51 cases with $m\ge3$, integral $N$ and $2\le a<a+1<N$. The largest $m$ is 56,917,632. Exact factorial valuations exhibit a common endpoint prime in each case. This proves $q\ge13$ without a weight cutoff. [Cremona integral-point tables at the revision used](https://github.com/JohnCremona/ecdata/tree/6fedd3ab463a70e52061ebc0b04832ad19fd2e01/intpts) · [LMFDB verification of integral points in Magma and SageMath](https://www.lmfdb.org/EllipticCurve/Q/Reliability)

In (5) the input point is integral. If $2\mid d$ and $8\mid q$, use $x'=x/4,y'=y/8$. Here $m$ is even since $q\mid m(m-2)$, and $t$ is even, so integrality is preserved. In all other cases use the original coordinates. Translate the resulting coordinates by integers to put the coefficient of $x^2$ in $\{-1,0,1\}$ and the linear coefficient of $y$ in $\{0,1\}$. Only exact matches to the published equations are used. Each tabulated integral $x$-coordinate gives at most one $m$, and the positive $a$ follows from its quadratic equation.

For example, at $q=1$ the transformation gives the curve $y^2+y=x^3+x^2$ with $x=m-1$. Its integral $x$-coordinates are $-1,0,1,2,21$. Only $m=22$, $N=441$, $a=98$ fits the cut conditions, and 2 divides both coefficients. [LMFDB curve 43.a1](https://www.lmfdb.org/EllipticCurve/Q/43/a/1)

## Enumeration through one million

The first program loops over $3\le m\le10^6$, tries every divisor of $m(m-2)$ as $N-1$ and tests whether $1+4mN$ is a square. There are 1,813 integral cut pairs, the largest at weight 9,896,020,401. Only six pass the parity test.

| $m$ | $N$ | $a$ | Least common endpoint prime |
| ---: | ---: | ---: | ---: |
| 84 | 493 | 203 | 3 |
| 126 | 187 | 153 | 3 |
| 3,712 | 7,421 | 5,248 | 5 |
| 120,872 | 187,561 | 150,568 | 3 |
| 422,010 | 654,841 | 525,689 | 3 |
| 588,800 | 7,476,801 | 2,098,175 | 3 |

Independently, loop over $1\le q\le10^6-3$, solve $(m-1)^2\equiv1\bmod q$ at its prime powers and combine by CRT. The bound $q<m-2$ follows from $N-1>m$. Set $N=1+m(m-2)/q$ and test the triangular-number discriminant. Both programs agree on the count, maximum and complete list of parity survivors. The valuation check factors $N$ and $N-1$ separately.

~~~sh
mkdir -p .work
c++ -std=c++17 -O3 audit_consecutive_cuts.cpp -o .work/audit_consecutive_cuts
.work/audit_consecutive_cuts 1000000 > .work/consecutive_cpp.log
python3 check_consecutive_cuts.py 1000000 --audit-log .work/consecutive_cpp.log
python3 check_consecutive_curves.py --download
lean ConsecutiveCutModels.lean
~~~

The Python programs use python-flint. The curve checker downloads the pinned 50 integral-point files into the ignored project directory and reuses them on subsequent runs.

Lean 4.36.0 checks the elimination of the weight, the short equation, its discriminant identity and the square-gap identity. [Lean source](ConsecutiveCutModels.lean)
