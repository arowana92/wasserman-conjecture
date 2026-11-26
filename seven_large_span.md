# Seven digit units at arbitrary positions for predecessor primes at least 11

Let $p\ge11$ be a prime divisor of $N-1$ with $s_p(N)=7$. Write the six positive digit positions, including multiplicity, as $a+h_i$, where $a\ge1$ and $0=h_1\le\cdots\le h_6=H$.

**If $H\ge10^8$, every three positive trinomial coefficients of weight $N$ have a common prime divisor.** At $p\ge29$ the [base-27 theorem](seven_base27.md) proves this for every $H$. Consequently only finitely many choices of $(p,h_2,\ldots,h_6)$ with $p\ge11$ can fail at any outer shift. The six positions here are arbitrary. The proof does not require them to come from a product or an arithmetic progression.

The [small-prime extension](seven_small_prime_span.md) also covers 5 and 7 for $H\ge10^{16}$.

Assume $11\le p\le29$. Use the [seven-unit reduction](seven_span_and_gaps.md#which-allocations-need-checking) at $B=p$. Write the binomial as $x=X^i F$ with $F=1+X^d$, the trinomial as $y=X^j(1+X^e+X^f)$, $0\le e\le f$, and write $z=X^s$. The assigned divisors require
$$U\mid F(p),y(p),\quad A\mid F(p),(y+z)(p),\quad C\mid y(p),(x+z)(p),\qquad 27p^{2H}\le U^3 A^2 C^2. \tag{1}$$
Pure-power block cases are already excluded in the reduction.

## Forcing a large cyclotomic factor

The arithmetic estimates for $AC$ in the [mixed cutoff proof](seven_base31.md#a-sparse-prime-of-the-weight) apply to these allocations without any second sparse prime. That extra hypothesis was used there to require $m=2AC$. For $H\ge10^5$ they give the stronger uniform estimate
$$AC<p^{197H/200}.$$
To see this, write $k=\deg\gcd(F,y+z)$. For $k=0,1$ use $3\sqrt{10}$ or $4\sqrt6(p+1)^{1/H}$ for $(AC)^{1/H}$, as in the mixed proof. For $1<k\le d/3$ and reduced norm $L\le2$, the common polynomial has value at most $4p^k$: there are at most two matching binomials $X^{j_i}+1$, and their monic lcm is their product divided by their common binomial factor when both occur. Consequently $AC\le4[3\cdot2^{2/3}p^{1/3}]^H$. At $L>2$ all four terms are congruent modulo $k$, giving $AC\le2[3\cdot6^{1/4}p^{1/3}]^H$, which is smaller. Here $C\le3^H$ since the alternative requires $k\le1$.

For $k>d/3$ use $k=d/3+2g$ and $A\le2^{(d-k)/g}p^k$. At $g\ge10$ this gives $(AC)^{1/H}\le3\cdot2^{1/10}p^{7/15}$. At $g\le9$ use $AC\le p^{18}[3\cdot2^{2/3}p^{1/3}]^H$. All these bounds are strictly below $p^{197H/200}$, checked at $p=11$, $H=10^5$, with decreasing quotient for larger $p,H$.

Thus $U>p^{H/100}$. The [binomial–trinomial gcd theorem](binomial_trinomial_gcd.md) bounds $U$ by $V p^{4H^{2/3}}$, where $V=1$ or $p^{2t}-p^t+1$ for $t=\gcd(d,e,f)$. For $H\ge10^8$ we have $4H^{2/3}\le 4H/460$. Since $1/100-4/460>1/800$, we need $V\ne1$ and $t>H/1600$. In particular $d/t\equiv3\bmod6$, and the other two gap quotients are 2 and 4 modulo six in some order.

The [seven-unit root classification](seven_large_bases.md#the-size-inequality) then gives $\gcd(F,y+z)=1$. Directly, a zero of $F$ allowing the four unit terms of $y+z$ to cancel would pair two of the terms of $y$ oppositely. Their differences are divisible by $2t$, whereas $v_2(d)=v_2(t)$. Also $C\le3^H$ since the three exponents of $y$ are distinct, by the [companion resultant bound](seven_base31.md#two-bounds-on-the-companion-factor).

If $d=3t$ and $H<6t$, the [short-span part of the cyclotomic argument](seven_refined_cutoffs.md#the-cutoff-at-64) gives $C=1$, $U<3p^{2t}$, $A\le3^t+1$ and $H\ge4t$. Since $3^t+1<p^t$, (1) already fails. In every other case $6t\le H$. It remains to control $A$ when the cyclotomic contribution to $U$ is large.

## A second gcd in two exponent variables

Put $D=H/t$, $v=s-j$, $n=d/t$, and $Y(S)=1+S^{e/t}+S^{f/t}$. Then $A$ divides the evaluations at $(S,R)=(p^t,p^v)$ of
$$1+S^n,\qquad Y(S)+R.$$
For **every** primitive integer direction $(u,w)$, set $S=p^t T^u,R=p^v T^w$. The two resulting Laurent polynomials have no common torus root. If $u=0$ the first is a nonzero constant. Otherwise a common root would have $S$ a root of unity and $T=\xi p^\lambda$, with $\xi$ a root of unity and $\lambda=-t/u$. This gives
$$p^{2(v+\lambda w)}=|Y(S)|^2=Y(S)Y(\overline S).$$
The right side is a nonzero algebraic integer with all conjugates real in $[0,9]$. As in the [rational-power argument](binomial_trinomial_gcd.md#locating-every-common-root), a rational power of a prime with these properties is an integer, and $p\ge11$ forces it to be one. Now $(t,v)$ is proportional to $(u,w)$. Its scalar is an integer by primitivity and nonzero since $t>0$. A common root along this homogeneous line would lift over $\mathbb C$ to a common root of $1+X^d$ and $Y(X^t)+X^v$: if $(t,v)=b(u,w)$ take $X^b=p^b T$ at the putative root. This contradicts $\gcd(F,y+z)=1$ after clearing monomials.

For each integer $Q\ge2$ choose $1\le q\le Q$ and an integer $r$ with $|qv/t-r|\le1/Q$. This is the elementary pigeonhole approximation of the fractional parts of multiples of $v/t$. In direction $(q,r)$, approximate the point $(t,v)$ by $(t/q)(q,r)$. The first coordinate is exact. The second has error at most $t/(qQ)$. The same absence of common roots holds for a nonprimitive direction, by passing to the corresponding power of $T$.

After stripping powers of $T$ the two degrees are $nq\le Dq$ and at most $2Dq+1$ respectively, since $|r|\le Dq+1$. Their nonzero resultant is divisible by $A$ away from $p$, by evaluating its identity at $T=1$. Every determinant monomial has the same $(q,r)$-weight, by the [weighted determinant argument](translated_triple_weights.md#a-weighted-determinant-and-an-explicit-bound). The absolute multiplicity sum is at most $2^{2Dq+1}4^{Dq}\le p^{2Dq}$, since $p\ge11$ and $Dq\ge1$. The only term with any residual error is $R$ in $Y(S)+R$. At most $nq$ columns use that polynomial, so the width of the powers of $p$ in the resultant is at most $nq\cdot t/(qQ)=d/Q$. Removing the lowest power of $p$ gives
$$ \log_p A\le 2D Q+H/Q. $$
The first polynomial has positive degree; the second also has positive degree because the distinct exponents in $Y$ are multiplied by $q>0$. Taking $Q=\lceil\sqrt t\rceil\le3\sqrt t/2$ gives $\log_p A\le4H/\sqrt t$.

## The size contradiction

Because $9^{16}<11^{15}$, we have $C^2\le p^{15H/16}$. Together with the bounds for $U,A$, this gives
$$\log_p(U^3 A^2 C^2)\le6t+12H^{2/3}+8H/\sqrt t+15H/16.$$
Here $6t\le H$. Since $t>H/1600\ge62500$ and $H^{1/3}\ge460$, the two error terms divided by $H$ sum to less than $12/460+8/250<1/16$. This contradicts (1).

Run the following with SymPy and python-flint installed. It checks the rational inequalities and 6,750 exact polynomial gcds on signed affine exponent lines, including constant-coordinate and homogeneous cases.

~~~sh
python3 check_seven_large_span.py
~~~
