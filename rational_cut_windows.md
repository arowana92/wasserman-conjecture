# Two binomial cuts near arbitrary rational proportions

Fix reduced fractions $0<r/s,u/t<1$ and $C>0$. Take any interior cuts $a,b$ at the same weight $N$ and set
$$x=sa-rN,\qquad y=tb-uN,\qquad M=N-1,\qquad R=r(s-r),\quad T=u(t-u).$$
Consider pairs with $|x|,|y|\le C\sqrt N$ whose two binomial coefficients share no prime of $NM$. Every cut of a positive trinomial row gives a binomial divisor, so the results also restrict the two nonendpoint rows in any counterexample.

* Except for the exact forms below, each interval of length $\sqrt X$ in $[X,2X]$ contains $O_{r,s,u,t,C}(1)$ possible weights. Through $X$ their total number is $O_{r,s,u,t,C}(X^{(\log_2\varphi)/2})$, where $\varphi=(1+\sqrt5)/2$ and the exponent is less than $0.348$.
* The exact forms have $a=rN/s$ and $b=uM/t$ or $uM/t+1$, or the interchanged description. In the indicated order they require $\gcd(s,t)=1$ and respectively $s\mid u$ or $s\mid t-u$. In particular when the fractions have the same denominator, the counting theorem applies to **all** the cut pairs in the stated windows. Across all the exact forms there are at most $O_{r,s,u,t}(X^{\log_2\varphi})$ possible weights through $X$.
* There is also an all-weight test. If
$$|xy|<N,\qquad R|(y+u)(y+u-t)|<M,\qquad T|(x+r)(x+r-s)|<M,$$
then only the stated exact forms can have endpoint gcd one. For example for two cuts near the *same* $r/s$, it suffices for a shared endpoint prime that both offsets have magnitude at most $c\sqrt N$ for any fixed $c<1/\sqrt R$ and all sufficiently large $N$. At the center a stronger [explicit interval](central_square_windows.md) holds at every weight.

## Local conditions and an integral conic

If the endpoint gcd is one, assign each full prime power of $N$ or $M$ to a cut avoiding its base. At a weight power the assigned cut is zero modulo the power, and at a predecessor power it is zero or one. Therefore $N\mid ab$ and $k=xy/N$ is an integer. Write $A=2r-s$, $B=2u-t$. At every predecessor power at least one of $(x+r)(x+r-s)$ and $(y+u)(y+u-t)$ vanishes. Thus
$$M\mid (x^2+Ax-R)(y^2+By-T).$$
Using $xy\equiv k\bmod M$ gives an integer $h$ with
$$hM=-T x^2-Ry^2+(Bk-AT)x+(Ak-RB)y+k^2+ABk+RT. \tag{1}$$
Under the window bounds $k,h$ range over finite sets depending only on the stated parameters.

If $x=0$, substitution in the divisibility before (1) gives $M\mid R(y+u)(y+u-t)$. The quotient is bounded under the window bounds. If it is nonzero it fixes $N$ as a quadratic polynomial of $y$ with nonzero leading coefficient. If it is zero we obtain the two exact values of $b$. Since $s\mid N$ and $t\mid M$, the denominators must be coprime. In the first exact form each $p^e\parallel N$ with $p\mid s$ fails to divide $a$ and hence divides $b=uM/t$; in particular $s\mid u$. In the second use the complementary cut $N-b=(t-u)M/t$. The case $y=0$ is symmetric. This also proves the explicit all-weight test, since $|xy|<N$ forces $k=0$.

If $k\ne0$, substitute $N=xy/k$ into (1). For each of the finitely many pairs $k,h$ this is a conic
$$T x^2+(h/k)xy+R y^2+(AT-Bk)x+(RB-Ak)y=k^2+ABk+RT+h. \tag{2}$$
In a dyadic interval $[X,2X]$, the absolute values of $x,y$ are both of order $\sqrt X$ since their product is $kN$. Solving (2) for $y$ gives at most two branches consisting of an affine function of $x$ plus or minus a constant times the square root of a polynomial of degree at most two. On each unbounded real branch as $x\to+\infty$ or $-\infty$, $y/x$ and $dy/dx$ approach the same number $z$ with $Rz^2+(h/k)z+T=0$. This follows directly by differentiating the quadratic formula: a square root of a linear polynomial contributes $O(\sqrt{|x|})$ and derivative $O(1/\sqrt{|x|})$, a square root of a quadratic with positive leading coefficient is linear to leading order, and an identically zero discriminant gives a line. Since $T>0$, $z\ne0$. The derivative of $N=xy/k$ consequently has absolute value $\gg\sqrt X$ on each relevant tail. An interval of weights of length $\sqrt X$ contains only $O(1)$ integer choices for $x$, hence only $O(1)$ pairs on these branches. The nonconstant quadratic in the $k=0$ case has the same spacing. Bounded branches and the finite initial portions of the tails are absorbed into the constant.

## Counting the eligible intervals

Partition a dyadic range into intervals $[Hn,H(n+1))$ with $H$ a power of two between $\sqrt X$ and $2\sqrt X$. At least one of the two cuts must be odd in coefficient. Suppose it is $a$. Its two sides have disjoint binary ones and their upper blocks $\lfloor a/H\rfloor,\lfloor(N-a)/H\rfloor$ add without carry to $n$. Write $n=sv+j$, $0\le j<s$. The first upper block is $rv+e$ and the other is $(s-r)v+j-e$, with $e$ in a fixed finite set because $|sa-rN|\le C\sqrt{2X}$. For fixed $e,j$ these are affine functions of one integer with positive slopes. A fixed shift of the parameter makes both offsets nonnegative and accounts for only finitely many initial inputs. The binary affine counting theorem bounds the number of accepted inputs up to $O(\sqrt X)$ by $O(X^{(\log_2\varphi)/2})$. Use the same argument if the other cut is odd. Combined with the spacing above and a sum over dyadic ranges, this proves the first count. [Binary affine counting theorem](binary_affine_counts.md)

For each exact form, restrict $N$ to a class modulo $st$. The two sides of either cut are now affine functions with positive integer slopes in a parameter of size $O(X)$. Apply that same binary bound directly to whichever cut is odd to obtain the exact-form count.

## Two complete trinomial rows

Fix two positive rational vectors $\rho,\sigma$ of length three, each summing to one. Suppose the two nonendpoint rows of a counterexample have
$$|A_i-\rho_i N|,\ |B_i-\sigma_i N|\le C\sqrt N\qquad(1\le i\le3).$$
Outside the two exact forms
$$A=N\rho,\quad B=M\sigma+e_j,\qquad\text{or}\qquad B=N\sigma,\quad A=M\rho+e_j,$$
there are at most $O_{\rho,\sigma,C}(X^{(\log_2\lambda_3)/2})$ such weights through $X$, where $\lambda_3^3=\lambda_3^2+1$ and the exponent is less than $0.276$. If the least common denominators of the two vectors have a common prime, the bound includes every such weight. The exact forms themselves contribute $O_{\rho,\sigma}(X^{\log_2\lambda_3})$ weights. This complements the [rigidity for deviations smaller than the square-root scale](projected_rows.md).

To prove this, apply the conic argument to all nine pairs of single parts. If any pair is outside the exact mixed forms, it gives $O(1)$ weights per interval of length $\sqrt X$. If all nine pairs have an exact mixed form, each has exactly one part at its exact weight proportion, since $0<\rho_i,\sigma_j<1$. This forces all three entries of one row to be exact and all three of the other to be a predecessor proportion plus zero or one. The sum of the latter corrections is one. The common denominators then divide $N$ and $M$ in some order.

For the count of eligible intervals, use the row that avoids 2. With $s$ the common denominator of its prescribed vector and $n=sv+j$ the upper block of the weight, the three upper parts are three affine functions of $v$ with fixed positive slopes and finitely many possible integer offsets. Their binary digits add without carry. Apply the three-part version of the binary affine count to $v=O(\sqrt X)$, or directly to an affine parameter of size $O(X)$ on the exact forms.

There is also a restriction on a row that need not be close to balanced. If $a$ is the middle part of a sorted row, the difference between its two successive gaps is $3a-N$ up to sign. Taking both rational fractions to be $1/3$ in the all-weight test shows that in a counterexample the two rows cannot both have successive gaps differing by at most $c\sqrt N$ for fixed $c<1/\sqrt2$ and all sufficiently large $N$. The lengths of the gaps themselves are unrestricted.

## Reproduction

Run `python3 check_rational_cut_windows.py`. The exact binomial calculation checks 334,958 retained-power divisibilities through weight 80, including 247,816 with an intermediate endpoint gcd. For every ordered endpoint-coprime cut pair through 400 it also tests all pairs of reduced fractions of denominators 2, 3, 4, 5, 7 and 9. This gives 3,360 nonzero-product conics, 156 nonconstant quadratic cases and 12 exact forms. One exact example at $N=46$ has $a=23$ and $b=10=2(N-1)/9$.
