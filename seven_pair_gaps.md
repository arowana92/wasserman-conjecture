# A binomial–trinomial gcd bound and the gaps of two digit pairs

For any odd prime $p$, any $n\ge1$ and $r,s\ge0$, define
$$G=\gcd(p^n+1,1+p^r+p^s),\qquad
D=\operatorname{odd}\big(\gcd(p^n+1,p^{|r-s|}+1)\big),$$
where $\operatorname{odd}$ removes the entire power of two. Then
$$G^3D^2<27p^{2n}. \tag{1}$$

**Application to seven units.** Suppose $p$ is odd, $p\mid N-1$, $s_p(N)=7$, and there is a counterexample at $N$. A row avoiding $p$ splits the units into two pairs and one triple, as in the [seven-unit reduction](seven_digit_units.md#two-pairs-at-an-odd-predecessor-prime). Suppose the unit at exponent zero belongs to the triple. Let $d,e$ be the gaps inside the two pairs, $g=\gcd(d,e)$, and $H$ the difference between the greatest and least of the six positive exponents. A necessary condition is
$$H< d+e-g. \tag{2}$$
In particular the two pairs cannot have equal gaps, even when the lower exponents differ. This assertion allows all six positive exponents to be distinct.

Indeed both $d/g,e/g$ are odd and the assigned weight scale $U$ divides $F=p^g+1$. The triple is $1+p^r+p^s$, so $U\mid G$ at $n=g$. Let $m=(N-1)/p^{v_p(N-1)}$ and take the predecessor factor $S$ assigned to this row after removing $p$. It is odd, since the companion must avoid 2. Each full power of $S$ divides one of the unshifted pairs, so $S\mid\operatorname{lcm}(p^d+1,p^e+1)$. At $\gcd(S,F)$ both pairs vanish, hence $p^r+p^s$ vanishes there by summing the three parts modulo that divisor of $N-1$. Thus $\gcd(S,F)\mid D$. Since the two binomials have gcd $F$,
$$S\le D\frac{p^d+1}{F}\frac{p^e+1}{F}.$$
The complementary size bound gives $27m^2\le U^3S^2$. By (1) and $(p^d+1)/F\le p^{d-g}$ we obtain $p^H\le m<p^{d+e-g}$, proving (2).

**If the unit at zero belongs to a pair**, let that pair be $1+p^u$, let the other pair have gap $d$, and put $g=\gcd(u,d)$. Then
$$H<u+g/3. \tag{4}$$
Indeed the third part divided by its lowest power of $p$ is $1+p^r+p^s$. The scale $U$ divides this and $p^g+1$, so (1) gives $U<3p^{2g/3}$. At every power in $S$ the first pair cannot be the exceptional part of residue 1, since that would make the power divide $p^u$. Thus $S\mid p^u+1$ and, by coprimality, $US\mid p^u+1$. Consequently $m\le (p^u+1)\sqrt U/\sqrt{27}<(p^u+1)p^{g/3}/3<p^{u+g/3}$. Together, (2) and (4) give tests using the exponent positions for every location of the unit digit.

## Congruent exponent positions

Put $F(X)=X^n+1$, $E(X)=1+X^r+X^s$. If $n\mid r-s$ with odd quotient, $E\equiv1\bmod F$, so $G=1$ and $D\le p^n+1$. If the quotient is even, $E\equiv1+2X^r\bmod F$ and $D=1$. From $2p^r=-1\bmod G$ and $p^n=-1\bmod G$ one gets $G\le2^n+1$ by taking $n$-th powers. Both cases give (1) for $p\ge3$.

Suppose $n\nmid r-s$ but $n$ divides one of $r,s$. Modulo $F$, $E$ becomes a monomial or $2+X^w$, where $g=\gcd(n,w)=\gcd(n,r-s)$. The former gives $G=1$. In the latter case $p^w=-2\bmod G$, whence $G\le2^\alpha+1$ at $\alpha=n/g\ge2$. Also $D\le B+1$ at $B=p^g$: a common divisor of two plus binomials above two forces both exponents divided by their gcd to be odd, and the gcd is then the smaller binomial $p^g+1$. Therefore
$$G^3D^2\le(2^\alpha+1)^3(B+1)^2\le (125/36)8^\alpha B^2<27B^{2\alpha}.$$
The last comparison holds already at $\alpha=2$, $B=3$, and the ratio decreases with either parameter.

## Three distinct positions modulo $n$

Let $h=\gcd(F,E)$ be monic over $\mathbb Q[X]$, and $k=\deg h$. Then
$$G\le 3^{(n-k)/2}|h(p)|. \tag{3}$$
Multiplication by $E$ on $\mathbb Q[X]/(F)$ has an integral matrix in the power basis. Since the positions are distinct modulo $n$, its columns each have three entries $\pm1$ and all others zero. Its rank is $n-k$ since its kernel consists of multiples of $F/h$. Its image contains the residue of the integral polynomial $h$, by Bézout and Gauss's lemma. Choose a maximal nonsingular minor, of absolute determinant $A$. Its columns in the full matrix form a basis of the image. Cramer's rule on the chosen rows shows that $Ah$ is an integral combination of those columns in the quotient. Evaluation at $p$ gives $G\mid Ah(p)$ since $F$ is monic. Hadamard's determinant bound proves (3).

If $k>0$, set $t=\gcd(n,r,s)$. At a common root $\zeta$, the three unit complex numbers $1,\zeta^r,\zeta^s$ sum to zero, so the last two are the distinct primitive cube roots. The order of $y=\zeta^t$ divides $\gcd(2n/t,3r/t,3s/t)$, a divisor of 6 by Bézout. Its order is both even and divisible by 3, so it is 6. Thus $n/t\equiv3\bmod6$ and $r/t,s/t$ are congruent to 2 and 4 in some order modulo 6. Conversely all $\zeta$ with $\zeta^t$ primitive sixth satisfy both equations under these congruences. The roots of $F$ are simple, giving $h=X^{2t}-X^t+1$. Here the 2-adic valuations of $n$ and $|r-s|$ differ, so $D=1$. If $n=3t$, then $F=(X^t+1)h$ and $E(p)\equiv3\bmod(p^t+1)$ since $r/t,s/t$ are even. Hence $G\le3h(p)<3p^{2t}$ and (1) follows. Otherwise $n\ge9t$ and (3) gives $G^3\le p^{3n/2+3t}\le p^{2n}$.

Finally suppose $k=0$, so $G\le3^{n/2}$. If $D=1$ this suffices. Otherwise, with $g=\gcd(n,r-s)$, the quotients $n/g,|r-s|/g$ are odd; the first is at least three since the positions are distinct modulo $n$. We have $D\le p^g+1$. If $n/g\ge5$, then $G^3D^2\le4p^{3n/2+2n/5}<27p^{2n}$.

When $n=3g$, $E$ takes the value 1 on every root of $X^g+1$, since $r-s$ is an odd multiple of $g$. On all $n$ roots of $F$ the sum of $|E|^2$ is $3n$ by orthogonality, using the three distinct positions. On the remaining $2g$ roots the sum is therefore $8g$. Their product of absolute values is at most $4^g$ by the arithmetic–geometric mean inequality. This product equals $|\operatorname{Res}(F,E)|$, a nonzero integer divisible by $G$. Hence $G^3D^2\le4\cdot64^g p^{2g}<27p^{6g}$ since $64<3^4$. This handles base three as well.

The [seven-unit verifier](check_seven_digit_units.py) compares the polynomial gcd with the root classification and checks the numerical bounds, including $r=s$, negative exponent differences, and both collision and distinct-position cases.
