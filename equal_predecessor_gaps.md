# Equal gaps are necessary at six predecessor digit units

**Theorem.** Let $p$ be a prime divisor of $N-1$. If $s_p(N)=6$, a counterexample requires a row of the following form, up to permutation:
$$
(1+p^u,\ p^i(1+p^d),\ p^j(1+p^d)),\qquad d\mid u,\quad u/d\ \text{odd}.
$$
Here $i,j,d,u$ are positive. The five positive digit exponents of $N$, counted with their multiplicities, have span strictly less than $3d/2$. The two equal differences must use four separate digit units. At most five units already suffice to prove the conjecture. The equal-gap restriction also gives [criteria using two sparse primes](#two-sparse-primes), including one prime from the weight and one from its predecessor.

In particular, take any positive $m$ coprime to $p$ with exactly five base-$p$ digit units. If its five positions have no two disjoint pairs at the same positive distance greater than two thirds of their total span, **every weight $1+p^a m$, $a\ge1$, satisfies the full conjecture**. For instance, at every prime $p$ and all $a,t\ge1$, this proves the assertion at
$$
1+p^a(1+p^t+p^{4t}+p^{9t}+p^{11t}).
$$
The ten positive pair distances before scaling by $t$ are distinct. Among five-unit multipliers $m<p^L$ coprime to $p$, at most $O(L^3)$ can admit even one counterexample as $a$ varies. Indeed the lowest position is zero. Three of the other four positions determine the fourth in the equal-distance equation, if all four appear in it. If one of them is outside the equation, there is still a nontrivial linear relation among the four. There are only finitely many choices of the four units and their pairing. The number of such multipliers without the equal-distance restriction has order $L^4$.

## The two bin bounds

Use the [six-unit spacing theorem](six_digit_gaps.md), which gives a row
$$
(1+p^u,\ p^i(1+p^d),\ p^j(1+p^e))
$$
with $d,e\mid u$ and both quotients odd. Assume $d>e$ and set $g=\gcd(d,e)$, $h=d/g$, $k=e/g$, $G=p^g+1$. The coprime integers $h,k$ are odd and the spacing theorem gives $d<3e+g$. Write $m=(N-1)/p^{\min(u,i,j)}$ in this proof. It is at least $p^d$ because the span of the five positive exponents is at least $d$.

We recall the necessary assignment conditions from that theorem. The row takes a full factor $U\ge3$ of $N$, with $U\mid G$. The prime powers of $m$ assigned to it fall into two products $S_d,S_e$, with
$$
\begin{split}
m&\le (U/3)^{3/2} S_d S_e,\\
S_d&\mid (p^d+1)/G,\qquad S_d\mid p^j(1+p^e)-1,\\
S_e&\mid (p^e+1)/G,\qquad S_e\mid p^i(1+p^d)-1. \tag{1}
\end{split}
$$
For example the powers in $S_d$ divide the first two parts, so the third part is 1 modulo those powers. For odd $p$, $U$ is odd. The three parts have the same 2-adic valuation since each of $u,d,e$ is an odd multiple of $g$. Their lowest binary ones collide. Thus the row cannot take the prime 2, and $U\le G/2$.

The following bounds apply to $S_d$ and, after interchanging indices, to $S_e$.

## Root-of-unity products

Let $v$ divide $d,e,j$, and put $n=d/v$, $s=e/v$, $b=j/v$. At $X=p^v$ the polynomials $X^n+1$ and $X^b(1+X^s)-1$ vanish modulo $S_d$. Thus $S_d$ divides their integer resultant, by the resultant Bézout identity. At a root $\zeta$ of the first polynomial, $(\zeta^s)^h=-1$ since $k$ is odd. In particular $|1+\zeta^s|\ne1$, because equality would make $\zeta^s$ a primitive cube root. The resultant is nonzero.

If $n\mid s$ then $h=1$ and $S_d=1$. If $n\mid b$ or $n\mid b+s$, reduction at $X=p^v$ gives either $S_d=1$ or $p^e\equiv-2$ or its inverse modulo $S_d$. Since $(p^e)^h\equiv-1$, the last two cases give $S_d\mid2^h-1$. Otherwise $0,b,b+s$ are distinct modulo $n$. Orthogonality on the roots of $X^n=-1$ gives mean square 3 for $X^b(1+X^s)-1$ there. The geometric mean of its absolute values is at most $\sqrt3$. In particular,
$$
S_d\le 2^{d/v},\qquad
S_d\le 3^{d/(2v)}\ \text{in the distinct-residue case}. \tag{2}
$$

For $p\ge5$, take $v=1$. If $h\ge5$, then $k\le h-2$ and $5^{h-3/2}>2^{2h-2}$, first by $5^7>2^{16}$ at $h=5$ and then by monotonicity. If $h=3$, then $k=1$, $S_e=1$ and $5^{3/2}>2^3$. In either case (2) contradicts $p^{gh}\le m<p^{3g/2}S_dS_e$.

## Base three

For $g=1$ there is no odd $U\ge3$ dividing 4. For $2\le g\le4$, (2) with $v=1$ and the aligned bound $2^h\le3^{d/2}$ always give $S_d\le3^{d/2}$, and similarly $S_e\le3^{e/2}$. Put $z=3^g$. Then $((z+1)/6)^{3/2}<z$ since $9\le z\le81$ and $(10/9)^3\,81<216$. Hence (1) gives $m<z^{1+(h+k)/2}\le z^h$.

For $g\ge5$, write $\mathrm e$ for Euler's number. We claim
$$
\log S_d\le d/\mathrm e+(h-1)\log 2. \tag{3}
$$
If $\gcd(g,j)\ge2$ use it as $v$ in (2), since $(\log2)/2<1/\mathrm e$. If $\gcd(g,j)=1$, group the roots of $X^d=-1$ by $X^g=\eta$. Their $j$-th powers run over all roots of $Y^g=\eta^j$. The absolute product over this group is $|1-(1+\eta^k)^g\eta^j|$. It is 1 for $\eta=-1$ and at most $2\max(1,|1+\eta^k|)^g$ otherwise. For $h>1$, $\eta^k$ runs over the $h$ roots of $Z^h=-1$, with mean $|1+Z|^2=2$. The inequality $\log^+ r\le r^2/(2\mathrm e)$ proves (3). When $h=1$ the bound follows from $S_d=1$.

Apply (3) on both sides. The difference between $gh\log3$ and the resulting upper bound for $\log m$ from (1) is at least
$$
g[(h-3/2)\log3-(2h-2)/\mathrm e]-(2h-4)\log2
 +\tfrac32\log(6/(1+3^{-g})).
$$
Using $\log3>1.09$, $\log2<0.70$, $1/\mathrm e<3/8$ and $6/(1+3^{-g})>3$, this exceeds $(0.34h-0.885)g-1.4h+4.3>0$ for $h\ge3$, $g\ge5$. This excludes unequal gaps in base three.

## Base two with large common gap divisor

One of $b=1,2$ admits an integer $a\ge0$ with $r=bj-ag$ satisfying $|r|\le g/3$. Use the shortest circular gap among the residues $0,j,2j$ modulo $g$. At $X=2^g$, $S_d$ divides both $X^h+1$ and
$$2^{\max(r,0)}X^a(1+X^k)^b-2^{\max(-r,0)}.$$
The resultant is nonzero. At a root $\eta$ of $X^h+1$, $\eta^k$ has order twice an odd number. For $r>0$ a zero would equate an algebraic integer with $2^{-r}$. For $r=0$ it would require $|1+\eta^k|=1$. For $r<0$, the integer $|1+\eta^k|^b=2^{-r}$ can only be 2 with $b=2$, because $|1+\eta^k|<2$. This would make $\eta^k=\pm i$, impossible by its order. The product bound therefore gives $S_d\le 2^{h(g/3+3)}$ and likewise for $S_e$.

If $g\ge29$ and $h\ge5$, then $3g/2+(2h-2)(g/3+3)<gh$ since $g(2h-5)>36(h-1)$. For $h=3$ use $S_e=1$ and $3g/2+3(g/3+3)<3g$. Both contradict (1). For the rest of the proof $g\le28$.

## Approximating two exponent ratios in base two

We give the details of the effective bound needed for the finite check. In the next argument write $n$ for either gap, $s$ for the other gap, $j_0$ for the power outside $1+2^s$ in the matching divisibility in (1), and $S$ for the corresponding bin. Thus
$$
S\mid 2^n+1,\qquad S\mid 2^{j_0}(1+2^s)-1,\qquad v_2(s)=v_2(n).
$$
Put $K=10$. Dirichlet's box argument gives $1\le a\le K^2$ and nonnegative $b,c$ with
$$
|a s-b n|,\ |a j_0-c n|\le n/K.
$$
Indeed use the $K^2+1$ successive multiples of the two ratios modulo 1 and $K^2$ equal boxes. Write $n=a l+r$ with $0\le r<a$ and $x=2^l$. Define
$$
F(X)=2^rX^a+1,\qquad
H(X)=2^C(2^{j_0-cl}X^c+2^{j_0+s-(b+c)l}X^{b+c}-1)
$$
where $C$ clears the powers of two in the denominators. Both vanish at $x$ modulo the odd integer $S$. The roots of $F$ are $2^{-r/a}\zeta$ with $\zeta^a=-1$. Put $w=j_0-cn/a$, $v=s-bn/a$. At such a root the last polynomial is $2^C$ times
$$
2^w\zeta^c+2^{w+v}\zeta^{b+c}-1. \tag{4}
$$

**Common roots.** If (4) vanishes, let the two summands before the minus sign be $X_1,X_2$. They sum to 1. Among the rational exponents $0,w,w+v$ the least must occur twice. To see this, divide the equation by the corresponding power of 2 and reduce at a prime ideal over 2 in a number field containing the terms. Positive rational powers of 2 vanish there and roots of unity are units.

If $v=0$, their ratio $\xi=\zeta^b$ satisfies $\xi^{n/g}=-1$, since $bn=as$ and $s/g$ is odd. Every algebraic conjugate of $1+\xi$ has the same absolute value $2^{-w}$, by (4). Thus all conjugates of $\xi$ have the same real part and its order has Euler totient at most two. Its order is twice an odd number, leaving 2 and 6. These give $|1+\xi|=0,\sqrt3$, respectively, neither a rational power of 2. So $v\ne0$.

Now one of $w,w+v$ is zero and the other is positive. The corresponding summand of magnitude one is a root of unity $\xi$. All conjugates of $1-\xi$ have the same absolute value, a rational power of 2 greater than 1. Again the order of $\xi$ has totient at most two, hence belongs to $\{1,2,3,4,6\}$. Only orders 2 and 4 remain from the absolute values $0,2,\sqrt3,\sqrt2,1$. At order 2 the summands are $-1,2$. Their ratio has argument $\pi$, so $\zeta^b=-1$ and $v_2(b)=v_2(a)$. But $s=bn/a\pm1$ contradicts $v_2(s)=v_2(n)$.

At order 4 the summands are $i,1-i$ or $-i,1+i$, in either order. Thus
$$
s=(b/a)n\pm\tfrac12,\qquad v_2(b/a)=-2. \tag{5}
$$
Indeed their ratio has argument $3\pi/4$ or $-3\pi/4$, and $\zeta$ has order with 2-adic valuation $1+v_2(a)$. Only two values of $\zeta^b$ qualify. Each is taken by at most $\gcd(a,b)\le a/4$ roots of $\zeta^a=-1$. At most $a/2$ roots are common to $F,H$.

**Size bound.** Let $B(X)$ be their primitive polynomial gcd in $\mathbb Z[X]$, of degree $t$, choosing its sign so its leading coefficient is positive. Write $F=B D$. The factors $B,D$ have leading coefficients that are powers of 2, and $D$ and $H$ are coprime since $F$ is squarefree. The resultant identity shows $S\mid B(x)\operatorname{Res}(D,H)$.

The odd part of $|B(x)|$ is at most $(2^{n/a}+1)^t$. In fact $B(x)$ is a rational power of 2 times the product of $2^{n/a}-\zeta$ over its roots, with the power's denominator dividing $a$. Each factor in the product is an algebraic integer. Raising the product to the $a$-th power therefore gives a rational algebraic integer, hence an integer, with precisely the odd part of $|B(x)|^a$.

Similarly the odd part of the absolute resultant is at most $(3\cdot 2^{2n/(aK)})^{a-t}$. In (4) divide by $2^{\min(0,w,w+v)}$. The resulting values are algebraic integers of absolute value at most $3\cdot 2^{2n/(aK)}$, since the three exponents have span at most $|w|+|v|\le2n/(aK)$. The resultant is their product over the roots of $D$ times a rational power of 2 with denominator dividing $a$, so the same integer argument applies.

Using $\log_2 3<1.6$, we conclude
$$
\log_2 S \le
\begin{cases}
1.6K^2+2n/K, &t=0,\\
1.6K^2+n/2+n/K, &t>0.
\end{cases} \tag{6}
$$
For the second line use $t\le a/2$, $\log_2(2^{n/a}+1)\le n/a+1$ and the bound for the remaining $a-t$ roots.

## A finite test independent of the digit positions

If at most one bin has common roots, (6) and (1) give $d\le42+320+d/2+3d/10$, so $d\le1810$. We test all $d\le2000$. Otherwise (5) applies at both bins. There are rational slopes $R,R'$ of reduced denominator at most 100 such that
$$
R=(e\pm1/2)/d,\qquad R'=(d\pm1/2)/e,\qquad v_2(R)=v_2(R')=-2.
$$
In particular $v_2(d)=v_2(e)=1$, and $RR'\ne1$. Thus $1/10000\le|RR'-1|\le1/(2d)+1/(2e)+1/(4de)<3/d$, using $d<4e$. This gives $d<30000$.

For $2000<d<30000$ loop over its permitted residue $2\bmod4$. A denominator $A\le100$ for $R$ must divide $2d$ and satisfy $A\equiv4\bmod8$. Generate $e$ by $2e\equiv\pm1\bmod (2d/A)$, then check the denominator of $R'$. Independently, enumerate both reduced slopes and solve the two linear equations for $d,e$ at each choice of signs. Both methods give the same 98 pairs after imposing the gap hypotheses.

For each pair let $L_d=(2^d+1)/G$, $L_e=(2^e+1)/G$. A bin $S_d$ divides
$$
B_d=\gcd(L_d,(1+2^e)^{2d}-1),
$$
since $1+2^e$ is a power of 2 modulo $S_d$ and $2^{2d}\equiv1$. Interchange $d,e$ for $B_e$. The necessary inequality is $27\,2^{2d}\le G^3 B_d^2 B_e^2$. Among the 434,968 initial pairs and the 98 large pairs, only $(d,e)=(18,14)$ passes, with $G=5$, $B_d=481$, $B_e=3277$. At this pair scan $j\bmod 36$, $i\bmod28$ in (1). The two maxima for the actual bins are 481 and 113. These violate the necessary inequality.

A second calculation enumerates the initial pairs by the coprime odd ratios $h,k$. It uses $2^{d-e}-1$ for the order test at both bins. Modulo $S_d$ this integer is $2^{d-e}(1+2^e)$, and modulo $S_e$ it is $-(1+2^d)$. Both tests therefore follow from (1).

This exhausts the unequal gaps. At $d=e$, $G=p^d+1$ and (1) gives $S_d=S_e=1$. If the positive positions have span $H$, then $p^H\le m\le((p^d+1)/3)^{3/2}<p^{3d/2}$, proving the stated span bound. The remaining companion row must avoid every prime of $m$.

## Two sparse primes

The full conjecture holds under either of the following hypotheses.

1. Two distinct primes of $N-1$ each have $s_p(N)\le6$.
2. There are primes $p\mid N-1$ and $q\mid N$ with $s_p(N)\le6$ and $s_q(N)\le5$. The upper bound on $s_q(N)$ can be raised to six whenever $p=2$ or $q=2$.

For the first assertion we need only consider six units at each base. In a counterexample the row $A$ for $p$ has its two parts $p^i(1+p^d),p^j(1+p^d)$ coprime to every prime of $N-1$ other than $p$, since $1+p^d\mid N$. It must carry at every such prime because avoiding one requires two zero residues. The second base must therefore be avoided by the companion. Applying the same theorem there leaves no row to avoid a third prime of $N-1$. But weights with at most two prime factors of $N-1$ already satisfy the [full conjecture](all_weights_endpoint_reduction.md).

For the second assertion take the $p$-avoiding row $A$, with $s_p(N)=6$. Write $P=p^{v_p(N-1)}$, $m=(N-1)/P$, and assign $p$ to $A$. Its assigned weight factor is $U\le1+p^d\le1+m$, and its assigned predecessor factor is $P$. Thus $UP<2N$. The complementary row must cover the product $N(N-1)/(UP)>(N-1)/2$. If this row also avoided $q$ with $s_q(N)\le6$, the [single-coefficient inequality at a weight prime](six_units_weight_prime.md) would give an available product strictly less than $(N-1)/2$. Indeed the product of all full powers at primes it avoids divides $N(N-1)/\gcd(N(N-1),K_C)$. Thus $A$ must avoid $q$.

With $s_q(N)\le5$, one of the three parts of $A$ is a pure power of $q$. The last two parts are divisible by $p$, so only $1+p^u=q^b$ is possible. For odd $p$ this forces $q=2$, while the equal lowest binary positions of $A$ force a carry. This also handles $q=2$ with six units.

For $p=2$ the weight is odd and $q$ is odd, so even the bound $s_q(N)\le6$ implies $s_q(N)\le5$ by parity. If $1+2^u=q^b$ with $b>1$, any odd divisor of $b$ greater than one would give a nontrivial odd factor of $q^b-1$ by a geometric sum. Hence $b$ is even. The powers of two $q^{b/2}-1,q^{b/2}+1$ differ by two, giving $q=3$, $b=2$, $u=3$. Five distinct positive binary positions have span at least four, so the span bound and $d\mid u$ force $d=3$. The five consecutive positive positions must then be $1,2,3,4,5$: the two gap-three pairs use the first and fourth and the second and fifth, leaving the middle position $u=3$. This is the checked weight 63.

If $b=1$, $q=1+2^u$ divides every part of $A$, forcing $d=u$ because $d\mid u$ and $q\mid1+2^d$. Again $d\ge3$ by the span bound. Carry-free addition at $q$ gives $1+s_q(2^i)+s_q(2^j)\le5$. For an exponent at least $u$, its positive digit sum is zero modulo $q-1=2^u$, hence at least eight. Otherwise both powers are single digits, and $1+2^i+2^j\ge7$ since $i,j$ are distinct and positive. This contradiction proves the assertions.

## Two separated groups and consecutive weights

Suppose the six unit positions are $0,a,b,t,t+c,t+f$ with $0<a\le b\le H_0$, $0\le c\le f\le H_0$ and $3t\ge5H_0$. Then a counterexample requires $(c,f)=(a,b)$. This reduces the sufficient separation in the [earlier two-group criterion](six_digit_gaps.md#two-groups-of-three-digit-units) from $3H_0$ to $5H_0/3$.

To see this, pair the units as in the theorem, with the unit at zero matched to $u$ and two remaining pairs of equal gap $d$. If $u$ is in the lower group, one of these pairs crosses between groups and the other lies in the upper group. If $u=d$, the upper pair puts the largest exponent at least at $t+d$, and the smallest positive exponent is at most $d$. Thus $t<3d/2\le3H_0/2$. If $u\ge3d$, the crossing gives $t\le H_0+d\le4H_0/3$. Both are impossible.

If $u$ is in the upper group and neither pair crosses, the two smaller positive positions form a pair. Their smaller position is at most $H_0-d$, so the span bound requires $t<H_0+d/2\le3H_0/2$. Hence both pairs cross. If $u\ge3d$, the smallest positive position exceeds $u-3d/2\ge3d/2$. It follows that $d<2H_0/3$, while a crossing requires $t\le H_0+d$. This again contradicts the separation. Thus $u=d$, and the three upper positions are the three lower ones translated by $d$.

In particular, **if $N\equiv2\bmod4$ has six binary ones, the full conjecture holds at one or both of $N-1,N$**. Indeed an exception at $N$ must have the [even six-bit form](six_units_weight_prime.md) with $1<j<k<D/2$. At $N-1$ the lower positions are $0,j,k$ and the upper ones are $D+1,D+j,D+k$. Their widths are at most $k$ and $D+1>2k$. The upper offsets $j-1,k-1$ differ from the required $j,k$.

At a predecessor with $p\parallel N-1$, the span bound also forces $u=d$. The least positive exponent is 1, so $u<1+3d/2<3d$. At least one of $i,j$ equals 1: otherwise $u=1=d$ and $\max(i+d,j+d)\ge3$, contradicting the span bound. Thus every remaining weight has the more specific form $(1+p^d)(1+p+p^j)$ with $1\le j<1+d/2$.

The polynomial verifier checks the odd-part inequalities directly by integer Sylvester determinants on 4,120 polynomial pairs, including 31 with a nonconstant common factor. It also checks 1,691 separated patterns, allowing repeated digit units.

~~~sh
# Python with python-flint installed
python3 check_unequal_predecessor.py --compare-slopes
python3 check_unequal_predecessor_cross.py
python3 check_gap_resultants.py
~~~
