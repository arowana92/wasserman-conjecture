# Two independent exponents at every prime

**Theorem.** Let $p$ be any prime and let $a,j$ be any positive integers. Set $P=p^a$, $Z=p^j$ and
$$N=(P^2+2P+2)Z+2P(P+1).$$
Every three positive trinomial coefficients of weight $N$ have a common prime divisor.

The [binary proof](positive_root_weights.md#the-binary-proof-at-all-exponent-pairs) handles $p=2$. The argument below covers every odd prime, including the diagonal $a=j$. It uses the shared prime base to replace a divisibility bound at $p^a$ with one at $p^{|nj-ma|}$ for a suitable pair of small integers $(n,m)$. For example, near $j=2a$ the second argument can be as small as $p$ even when $P$ has arbitrarily many digits.

## The two divisors to bound

Write $B(X)=X^2+2X+2$, $C(X)=2X(X+1)$ and $M=N-1$. At any counterexample the [endpoint theorem](all_weights_endpoint_reduction.md) supplies $(1,1,N-2)$. One of the other rows $A$ avoids $p$, meaning that $p$ does not divide its coefficient. Split the full prime powers of $N=UV$ and $M=ST$ between the two nonendpoint rows, assigning each prime to a row which avoids it and assigning $p$ to $A$. The powers in $U,S$ go to $A$. Thus $U$ divides every $A_i$ and, at each power in $S$, two parts vanish and the third equals 1 modulo that power. The three parts of the companion row divided by $V$ sum to $U$. Their product is a multiple of $T^2$ and is at most $(U/3)^3$. Hence
$$27M^2\le U^3 S^2. \tag{1}$$

For odd $p$ the exact $p$-power in $N$ is $p^{\min(a,j)}$. Indeed the lowest coefficient is 2 when the exponents differ, and $N=P(P+2)^2$ when they agree. Put $U'=U/p^{\min(a,j)}$ and
$$h(a,j)=4a+2j-3\min(a,j).$$
Since $M>P^2Z$, (1) requires $27p^{h(a,j)}<(U')^3 S^2$.

By Kummer's criterion, $A_i=b_i(P)Z+c_i(P)$ for a coefficientwise distribution of $B,C$ among three positive parts. Overlapping digit positions are allowed. Merge identical base-$p$ units to form the digits, then give each original unit to the part receiving the merged unit. The [two-block minor calculation](positive_root_weights.md#two-independent-exponents-at-every-prime) supplies nonzero integer polynomials $R_0,R_1$ for each of the 575 distributions, with
$$U'\mid R_0(P),\qquad S\mid R_1(P).$$
Here $R_0$ is the bound with the previously restored factor $2X$ removed. To see why these same bounds work in both orders of the exponents, consider the three congruences $b_i(P)Z+c_i(P)-\delta_i=0$ modulo the assigned divisor, with $\delta=0$ for $U'$ or a unit vector for a bin of $S$. Eliminating $Z$ gives the two-by-two minors. Pure terms supply additional bounds, as $Z$ is invertible. Factors of $P$ can be canceled throughout. Taking a polynomial gcd and integer resultants bounds the common divisor, and the three pairwise coprime bin divisors together divide a polynomial lcm. This proves divisibility even at an evaluation where a bound vanishes.

## Projecting the exponent pair

Let $W=U'$ or $S$, and put $\sigma=0$ or 1, respectively. Then $p$ is invertible modulo $W$, $W\mid R(P)$ with $R=R_\sigma$, and $B(P)Z+C(P)\equiv\sigma\pmod W$. For a primitive direction $v=(n,m)$ with $n>0$, write $R=cR_*$ with $c$ its positive content, put $d=\deg R$, and form the integer polynomial
$$D_v(Y)=c\operatorname{Res}_X\!\left(R_*(X),\;Y X^m B(X)^n-(\sigma-C(X))^n\right).$$
For a constant $R$ one may simply take $D_v=c$. The resultant has degree at most $d$ in $Y$ by the Sylvester determinant. Put $e=nj-ma$. If $e\ge0$, the two equations modulo $W$ and the integral resultant identity give $W\mid D_v(p^e)$. Multiplying the identity by $c$ accounts for a nonprimitive $R$. If $e<0$, work modulo $W$ with $Y=(p^{-e})^{-1}$, then clear the denominator. In this case $W$ divides the evaluation of $Y^dD_v(1/Y)$ at $Y=p^{-e}$. The vertical direction $(0,1)$ uses $R$ at $p^a$ directly. The verifier checks every projected polynomial is nonzero.

Start with rays $(1,0),(1,1),(0,1)$. In each sector, subdivide at $v+w$ until each neighboring pair satisfies
$$h(v+w)\ge 3\deg R_0+2\deg R_1+2.$$
All neighboring determinants have absolute value one. The procedure terminates since $h(a,j)\ge a+j$ and the coordinate sum of the new mediant strictly increases along each branch. Interchanging neighbors when necessary, each integer exponent pair in the sector has a representation
$$(a,j)=(b+\ell)v+\ell w,\qquad b,\ell\ge0,\quad b,\ell\in\mathbb Z.$$
The projection argument at $v$ has absolute exponent $\ell$ by the determinant identity. Also $h(a,j)=b h(v)+\ell h(v+w)$ since each sector stays on one side of the diagonal.

Let $D_0,D_1$ be the oriented projections for this sector. Suppose neither vanishes at $p^\ell$. Write $E=3\deg D_0+2\deg D_1$ and $L=\lVert D_0\rVert_1^3\lVert D_1\rVert_1^2$, where the norm is the sum of absolute coefficients. The required inequality and the divisibilities imply
$$27p^{b h(v)+\ell(h(v+w)-E)}<L. \tag{2}$$
For $\ell=0$ one may replace $L$ by $|D_0(1)|^3|D_1(1)|^2$. Both coefficients in the exponent of (2) are positive. Taking $p\ge3$ therefore gives a finite list of pairs, with an integer upper bound on $p$ for each. If either projection vanishes, the complete positive integer root lists in the calculation have union $\{1,4,16\}$. Only 1 is possible for an odd prime, giving $\ell=0$. Interior rays occurring this way have $j/a$ equal to 1, 2 or 4.

## Checking the remaining cases

For a finite pair $(a,j)$ with prime ceiling at most 1,000, the verifier evaluates the row directly. It checks that (1) fails even with the numerical upper bounds
$$U\le\gcd_i A_i,\qquad S\le\gcd\!\left(M,\operatorname{lcm}_{i<k}\gcd(A_i,A_k)\right).$$
For larger ceilings and the three zero directions it checks the whole corresponding ray. Write $(a,j)=g(n,m)$ in lowest terms and use $t=p^g$. Then $A_i(t)=b_i(t^n)t^m+c_i(t^n)$ and $M(t)=\sum_i A_i(t)-1$. Removing powers of $t$ from the polynomials $A_i$ preserves divisibility by $U'$ and by the bins of $S$. Polynomial gcds and integer resultants give a bound $u_*(t)$ for $U'$ from the three parts and the weight. The same procedure on the predecessor and each pair of parts, followed by a polynomial lcm, gives a bound $s(t)$ for $S$. Put $u(t)=t^{\min(n,m)}u_*(t)$. The polynomial
$$729M(t)^4-u(t)^6s(t)^4$$
has positive constant term and nonnegative coefficients after a translation $t=t_0+z$ in every one of these cases. Thus it is positive for all $t\ge t_0$. Direct numerical gcds handle every integer $3\le t<t_0$. The largest translation used is 24.

The verifier checks 2,698 oriented cones across 575 shapes. It makes 18,675 direct prime and row evaluations, checks 147 shape and ray pairs symbolically, and makes 432 small integer evaluations on those rays. No evaluation satisfies the necessary bound. This completes the proof for all odd primes, and hence the theorem.

With SymPy and python-flint installed, run the following from the repository root.

~~~sh
python3 check_exponent_fan.py
python3 check_positive_root_weights.py
~~~
