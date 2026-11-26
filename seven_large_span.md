# Seven digit units at arbitrary positions for predecessor primes at least 11

Let $p\ge11$ be a prime divisor of $N-1$ with $s_p(N)=7$. Write the six positive digit positions, including multiplicity, as $a+h_i$, where $a\ge1$ and $0=h_1\le\cdots\le h_6=H$.

**If $H\ge10^{14}$, every three positive trinomial coefficients of weight $N$ have a common prime divisor.** At $p\ge29$ the [base-27 theorem](seven_base27.md) proves this for every $H$. Consequently only finitely many choices of $(p,h_2,\ldots,h_6)$ with $p\ge11$ can fail at any outer shift. The six positions here are arbitrary. The proof does not require them to come from a product or an arithmetic progression.

Assume $11\le p\le29$. Use the [seven-unit reduction](seven_span_and_gaps.md#which-allocations-need-checking) at $B=p$. Write the binomial as $x=X^i F$ with $F=1+X^d$, the trinomial as $y=X^j(1+X^e+X^f)$, $0\le e\le f$, and write $z=X^s$. The assigned divisors require
$$U\mid F(p),y(p),\quad A\mid F(p),(y+z)(p),\quad C\mid y(p),(x+z)(p),\qquad 27p^{2H}\le U^3 A^2 C^2. \tag{1}$$
Pure-power block cases are already excluded in the reduction.

## Forcing a large cyclotomic factor

The arithmetic estimates for $AC$ in the [mixed cutoff proof](seven_base31.md#a-sparse-prime-of-the-weight) apply to these allocations without any second sparse prime. That extra hypothesis was used there to require $m=2AC$. For $H\ge1000$ the estimates give, by themselves,
$$ AC\le 2^{-H/101}p^H\le p^{(1-1/600)H}. $$
Indeed all eight uniform bounds were checked after multiplication by $2^{1/101}$. The remaining two-matching case has $g=1$, $d<300$, $k=d/3+2$. Here $(AC)^{1/H}\le3\cdot2^{300/1000}p^{102/1000}<p/2^{1/101}$ at $p=11$ and the ratio decreases thereafter. The second inequality uses $p<2^5$.

Thus $U>p^{H/900}$. The [binomial–trinomial gcd theorem](binomial_trinomial_gcd.md) bounds $U$ by $V p^{13H^{2/3}}$, where $V=1$ or $p^{2t}-p^t+1$ for $t=\gcd(d,e,f)$. For $H\ge10^{14}$ we have $13H^{2/3}<H/1800$, so $V\ne1$ and $t>H/3600$. In particular $d/t\equiv3\bmod6$, and the other two gap quotients are 2 and 4 modulo six in some order.

The [seven-unit root classification](seven_large_bases.md#the-size-inequality) then gives $\gcd(F,y+z)=1$. Directly, a zero of $F$ allowing the four unit terms of $y+z$ to cancel would pair two of the terms of $y$ oppositely. Their differences are divisible by $2t$, whereas $v_2(d)=v_2(t)$. Also $C\le3^H$ since the three exponents of $y$ are distinct, by the [companion resultant bound](seven_base31.md#two-bounds-on-the-companion-factor).

If $d=3t$ and $H<6t$, the [short-span part of the cyclotomic argument](seven_refined_cutoffs.md#the-cutoff-at-64) gives $C=1$, $U<3p^{2t}$, $A\le3^t+1$ and $H\ge4t$. Since $3^t+1<p^t$, (1) already fails. In every other case $6t\le H$. It remains to control $A$ when the cyclotomic contribution to $U$ is large.

## A second gcd in two exponent variables

Put $D=\lfloor H/t\rfloor\le3600$, $v=s-j$, $n=d/t$, and $Y(S)=1+S^{e/t}+S^{f/t}$. Then $A$ divides the evaluations at $(S,R)=(p^t,p^v)$ of
$$1+S^n,\qquad Y(S)+R.$$
For **every** primitive integer direction $(u,w)$, set $S=p^t T^u,R=p^v T^w$. The two resulting Laurent polynomials have no common torus root. If $u=0$ the first is a nonzero constant. Otherwise a common root would have $S$ a root of unity and $T=\xi p^\lambda$, with $\xi$ a root of unity and $\lambda=-t/u$. This gives
$$p^{2(v+\lambda w)}=|Y(S)|^2=Y(S)Y(\overline S).$$
The right side is a nonzero algebraic integer with all conjugates real in $[0,9]$. As in the [rational-power argument](binomial_trinomial_gcd.md#locating-every-common-root), a rational power of a prime with these properties is an integer, and $p\ge11$ forces it to be one. Now $(t,v)$ is proportional to $(u,w)$. Its scalar is an integer by primitivity and nonzero since $t>0$. A common root along this homogeneous line would lift over $\mathbb C$ to a common root of $1+X^d$ and $Y(X^t)+X^v$: if $(t,v)=b(u,w)$ take $X^b=p^b T$ at the putative root. This contradicts $\gcd(F,y+z)=1$ after clearing monomials.

Use the elementary one-coordinate approximation to $(t,v)/K$ with $K=\max(t,|v|)\le H$. For each integer $Q\ge2$ it supplies $1\le q\le Q$, a primitive $(u,w)$ with norm at most $q$, and real $\alpha$ with $\|(t,v)-\alpha(u,w)\|_\infty\le H/(qQ)$. One coordinate is exact up to the scale and pigeonholing the other coordinate gives denominator at most $Q$.

Strip the lowest $T$ powers from the two polynomials. Their degrees are at most $Dq,2Dq$. Their nonzero resultant is divisible by $A$ away from $p$, since its integral Laurent identity can be evaluated at $T=1$. Every determinant monomial has the same $(u,w)$-weight, by the [weighted determinant argument](translated_triple_weights.md#a-weighted-determinant-and-an-explicit-bound). The multiplicity sum is at most $2^{2Dq}4^{Dq}\le p^{4Dq}$. The residual exponent contribution from a single entry lies between $-DH/(qQ)$ and $DH/(qQ)$. There are at most $3Dq$ entries in a product. Removing the lowest power of $p$ gives
$$ \log_p A\le4D Q+6D^2H/Q. $$
If either polynomial has degree zero, use that polynomial instead, with the same bounds. Taking $Q=\lceil\sqrt{DH}\rceil$ yields $\log_p A\le11D^{3/2}\sqrt H$ at the stated heights.

## The size contradiction

Because $9^{16}<11^{15}$, we have $C^2\le p^{15H/16}$. Together with the earlier bounds for $U,A$, this puts
$$ \log_p(U^3 A^2 C^2)\le 6t+39H^{2/3}+22D^{3/2}\sqrt H+15H/16. $$
If $t\le H/12$, then $6t\le H/2$ and $D\le3600$. The middle two terms divided by $H$ are at most $39/40000+22\cdot3600\cdot60/10^7<9/16$. If $t>H/12$, then $D\le12$ and the corresponding sum is at most $39/40000+22\cdot12\cdot4/10^7<1/16$. Here $6t\le H$ by the earlier case split. Both cases give $U^3A^2 C^2<p^{2H}$, contradicting (1).

Run the following with SymPy and python-flint installed. It checks the rational inequalities and 6,750 exact polynomial gcds on signed affine exponent lines, including constant-coordinate and homogeneous cases.

~~~sh
python3 check_seven_large_span.py
~~~
