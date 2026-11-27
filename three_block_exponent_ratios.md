# The only separated exponent ratio for three blocks

Fix a prime $p$ and positive integers $D,m,C$ with $C\ge2$. Let $X=\max(D,m,C)$, $n=\lceil\log_p(2X+2)\rceil$. Suppose
$$ N=D p^a+m p^b+C,\quad p\mid N(N-1),\quad
 a>16384(n+1),\quad b>4n,\quad a-b>4n. $$

**Theorem.** At a counterexample at least one of the following must hold.

1. In the row avoiding $p$, the three digit-block columns have the form $g z,h z,l z+e$, with $z$ a primitive positive integer triple, $g,h,l$ positive integers and $e$ zero or one coordinate unit. Each column adds without a base-$p$ carry. In particular $\sum z_i\ge3$ divides $D,m,C-\sum e_i$.
2. The exponents satisfy $|a-2b|\le6n$. There are two distinct targets among $0,e_1,e_2,e_3$ on the same affine plane of block columns, as described below.

In particular **coprime $D,m$ give a proof throughout the separated range outside the band $|a-2b|\le6n$**. The conclusion also holds outside this band whenever both $\gcd(D,m,C)$ and $\gcd(D,m,C-1)$ are below three.

**The entire separated range works**, including the band, if $D$ is one power of $p$, if the base-$p$ digit sum of $m$ is at most two, or if $C$ is one power of $p$. A power here may have exponent zero when consistent with the hypotheses. For instance this covers every $2^a+m2^b+C$ in the displayed range with completely arbitrary $m\ge1$ and $C\ge2$.

## The four congruences

Choose the row assigned $p$ in the [endpoint reduction](all_weights_endpoint_reduction.md), and split its parts into $d_i p^a+m_i p^b+c_i$. The three nonnegative columns sum to $D,m,C$ and separately split their digits without carries. Let $U,S_1,S_2,S_3$ be the products assigned to this row at the weight and at the three predecessor residue-one positions. The [complementary row](projected_rows.md) gives $27(N-1)^2\le U^3(S_1S_2S_3)^2$. Each actual part $y_i$ has $y_i^5>27(N-1)^2$, hence $y_i>X$. Tie the factors to the targets $e=0,e_1,e_2,e_3$, and let $W_e$ denote the part of the corresponding factor prime to $p$. All powers of $p$ removed together are at most $p^n$ since $2\le C<p^b$.

Put $B=2p^{128n+a/64}$ and call a factor $W=W_e$ large if $W>B$. For such a factor,
$$ W\mid d_i p^a +m_i p^b+c_i-e_i\quad(1\le i\le3), $$
where $e_i$ in the right-hand expression denotes coordinate $i$ of the chosen target. If the first two columns are dependent, a nonzero minor of the first and last adjusted columns would be divisible by $W$ yet at most $2X^2$. Thus all three columns are collinear. The first gives a primitive nonnegative direction $z$. A zero coordinate would make an actual part at most one. All multipliers are integral by primitivity, and $C-\sum e_i>0$ forces the last multiplier positive. This is alternative 1.

Otherwise the first two columns are independent. The determinant with the adjusted last column $c-e$ is divisible by $W$ and at most $6X^3$, so it vanishes. Cramer's rule gives a unique rational point $(x,y)=(I/\Delta,J/\Delta)$ with $d_i x+m_i y+c_i=e_i$ and
$$
 W\mid\Delta p^a-I,\quad W\mid\Delta p^b-J,\qquad
 0<|\Delta|\le 2X^2,\quad |I|,|J|\le2X^2.
$$
Neither numerator is zero. Every large target is therefore on the affine plane $c+\operatorname{span}(d,m_{\rm col})$, where $m_{\rm col}=(m_1,m_2,m_3)$. This plane contains at most two targets. Containing all three unit vectors would make its direction plane have coordinate sum zero, whereas the first column sums to $D>0$. Containing zero and two unit vectors would make both columns vanish on the remaining coordinate, giving an actual part at most $X$.

## Clearing denominators at the modulus

Use the same Dirichlet approximant $1\le u\le64$, $0\le v\le u$, $|va-ub|\le a/64$ for all targets. The two Cramer congruences imply $x^v=p^{va-ub}y^u$ as a rational identity. Indeed compare $I^v\Delta^{u-v}$ and $p^{va-ub}J^u$, clearing a negative power of $p$. The difference is divisible by $W$ and its absolute value is at most $B$. The endpoint approximants $v=0,u$ would give $b\le2n$ or $a-b\le4n$, respectively, from the $p$-valuations of $x,y$. Write $s/d$ for the reduced fraction $v/u$, so $0<s<d$. This denominator is shared by all large targets.

On absolute values, valuations outside $p$ give $|x|=p^j t^d$, $|y|=p^k t^s$ for a positive reduced fraction $t=A/H$ with numerator and denominator prime to $p$. Its numerator and denominator are at most $p^{2n}$, and $|j|,|k|\le2n$. The valuation at $p$ gives $(a-j,b-k)=g(d,s)$ for a positive integer $g$. These assertions follow prime by prime from $s\,v_q(|x|)=d\,v_q(|y|)$ at $q\ne p$ and from $s(a-j)=d(b-k)$.

Here the Cramer denominator can share primes with $W$. Let $h=\gcd(|\Delta|,|I|,|J|)$ and $W'=W/\gcd(W,h)$. Divide the triple $(\Delta,I,J)$ by $h$. The two congruences still hold modulo $W'$ for the resulting primitive triple. A prime dividing $W'$ and its new first entry would divide the other two entries as well, impossible. Hence $x,y$ have invertible denominators modulo $W'$ and reduce to the units $p^a,p^b$. Both $A$ and $H$ are invertible there too. For $Z=p^g/t$ in this unit group we have
$$ Z^d=\operatorname{sgn}(x),\qquad Z^s=\operatorname{sgn}(y). $$
Since $d,s$ are coprime, Bézout gives $Z=\pm1\bmod W'$. Therefore $W'$ divides one of $Hp^g\pm A$, both nonzero since $g>0$ and $p\nmid A$. In particular
$$ W\le |\Delta|(H p^g+A) \le 2p^{a/d+6n}. \tag{1} $$

Let $h_0$ be the sum of the exponents in $U^3(S_1S_2S_3)^2$ belonging to large factors. There are at most two such factors, so $h_0\le5$. With no large factor the product is at most $p^{3n}B^9<p^{2a}$. If $d\ge3$, the logarithm to base $p$ of the product is at most
$$ (5/3+9/64)a+1185n+14<2a $$
by (1), even overcounting the small-factor exponents by nine. If $d=2$ but there is only one large factor, replace $5/3$ by $3/2$ in the same upper bound. Each contradicts the complementary-row inequality. Thus $d=2$, $s=1$ and there are exactly two large targets. The valuation relation gives $a-2b=j-2k$, so $|a-2b|\le6n$.

## How sparse columns exclude the last ratio

For either of the two large targets $W'>2$. From $Z^2=\operatorname{sgn}(x)$ and $Z=\pm1\bmod W'$ we get $x>0$. Since $Dx+my+C=\sum e_i$ and $C\ge2$, we have $y<0$. Put $t_0=-y>0$ and $q=p^{a-2b}$. Then $x=q t_0^2$, since $|x|=p^j t^2$, $|y|=p^k t$ and $j-2k=a-2b$. The targets are therefore two distinct points on the curve
$$ e_i=d_i q t_0^2 -m_i t_0+c_i \quad (i=1,2,3). \tag{2}$$
Their parameters differ. A difference of two distinct vertices among zero and the coordinate units has at most one strictly positive and at most one strictly negative coordinate.

* If $D$ is one $p$-power, two $d_i$ vanish. At these positions $m_i>0$, since otherwise the actual part is at most $X$. Both coordinates in (2) strictly decrease with $t_0$, impossible at the two targets.
* If $m$ has digit sum at most two, some $m_i=0$. There $d_i>0$ by the same part bound. Since $x>0$, that target coordinate is positive at both points. The only permissible target with that coordinate positive is its unit vector.
* If $C$ is one $p$-power, two $c_i$ vanish. At both positions $d_i>0$, since otherwise (2) is negative or the actual part vanishes. Dividing the coordinate by $t_0$ gives the strictly increasing function $d_i q t_0-m_i$. It is nonnegative at both target parameters. Thus as $t_0$ increases, both original coordinates strictly increase, again impossible.

Each sparsity condition also rules out alternative 1: its corresponding digit column has three positive entries adding without a carry and therefore needs at least three digit units. This proves the stated exclusions.

## A quadratic example with two large divisors

Take $D=37,m=110,C=73$, $p=2$ and $a=2b$, $T=2^b$. The three digit columns are $(1,4,32),(2,12,96),(1,8,64)$. The corresponding quadratic part polynomials take values $(0,0,0)$ at $T=-1$ and $(1,0,0)$ at $T=-2$. Thus $T+1$ divides all three parts and $T+2$ divides their differences from the first coordinate unit. The columns split their bits freely, as do the evaluated parts for $b\ge7$. This verifies the two large congruences in alternative 2 at a carry-free row.

~~~sh
python3 check_three_block_ratios.py
~~~
