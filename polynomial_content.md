# Polynomial weights with a common coefficient factor

Let $H\in\mathbb Z[X]$ be primitive and irreducible over $\mathbb Q$, of degree $n\ge1$, with nonnegative coefficients and positive leading coefficient. Let $G=mH$ with integer $m\ge1$, let $P=p^a$ for a prime $p$ and $a\ge1$, and put $R=(1+\sum_i G_i^2)^{1/2}$. Assume $P>R^{15}$.

Call a positive integer $s$ excluded here if $s<364$, if it has a full prime power $w$ with $3w\ge s$, or if it is at most 1000 and is absent from the [sixteen possible shifted denominators](shifted_directions.md). Each condition independently proves that $s$ cannot be the denominator of a shifted row in a counterexample.

**Theorem for $PG(P)+1$.** Assume $XG+1$ is irreducible over $\mathbb Q$, and every divisor of $m$ is excluded here. The conjecture holds at $PG(P)+1$ if any of the following holds:

1. $m\le6$.
2. $H(0)\ne1$.
3. $m\ge7$ and $H_{i+1}>(m-4)H_i$ for some $0\le i<n$.
4. $m=r^e\ge7$ is a prime power and either $r\mid P$ or $jP\equiv1\pmod m$ for some $1\le j\le4$.

In condition 4 the divisor hypothesis is automatic, and in conditions 2 and 3 it is automatic for every $m\le363$ or $m=d r^e$ with $1\le d\le3$. Indeed every divisor of the latter either divides $d$ or has a full $r$-power of size at least a third of that divisor. These results allow arbitrary polynomial degree and arbitrary coefficient size. For fixed eligible $m,n$, conditions 2 and 3 leave only finitely many polynomials $H$ to consider, since the remaining coefficients satisfy $H_0=1$ and $H_i\le(m-4)^i$ when $m\ge7$.

**Theorem for $PG(P)$.** Assume $XG-1$ is irreducible. The conjecture holds at $PG(P)$ if $m/p^{v_p(m)}\le2$, or if some nonzero coefficient $G_i$ has base-$p$ digit sum at most two. In particular $m=1,2$ work for every prime base.

## Two ways a resultant can vanish

Use the [prime-power assignment and resultant estimates](polynomial_prime_powers.md#proof-by-specialization). In a counterexample choose a row $A$ avoiding $p$. Its three parts are $P L_i(P)+\delta_i$, where the nonnegative integer polynomials $L_i$ sum coefficientwise to $G$. At the shifted weight take $\delta=(0,0,1)$, at the unshifted weight take $\delta=0$. All $L_i$ are nonzero since the parts of $A$ must have a common divisor $U\ge3$. The coefficient splits at each power of $X$ are carry-free in base $p$, since $P>\max G_i$.

At the shifted weight $U\le R^{2n+1}$ by irreducibility of $XG+1$. Write $PB$ for the product of full powers of $PG(P)$ assigned to $A$. Group the powers left in $B$ by the slot $k$ whose actual part is 1 at that prime power. Within each group its product divides $G(P)$ and the evaluations of three polynomials. For $\delta_i=[i=k]$ take $L_i$ after cancelling $P$, and otherwise take $X L_i+\delta_i-[i=k]$. At the base $p$, the exceptional slot is 3 and the cancelled power $P$ was available in full. If one defining polynomial is coprime to $H$ over $\mathbb Q$, its resultant with $G$ bounds that group by $R^{2n+1}$. This follows from the Euclidean coefficient norms, both at most $R$, and degrees at most $n,n+1$.

If this works in all three groups then $B\le R^{6n+3}$, while the [assignment inequality](polynomial_prime_powers.md) gives
$$27G(P)^2\le U^3 B^2\le R^{18n+9}.$$
Since $G(P)\ge P^n$ and $(18n+9)/(2n)\le15$, this contradicts the hypothesis on $P$. Thus some group has all three defining polynomials divisible by $H$.

If $k=3$, we have $L_i=t_i H$, where the $t_i$ are positive integers of sum $m$. Their integrality uses primitivity and Gauss's lemma. The shifted denominator of $A$ is $m/\gcd(t_1,t_2,t_3)$, a divisor of $m$, excluded by hypothesis.

If $k\ne3$, rename the unshifted slots so $k=1$. We have $H\mid XL_1-1,L_2,XL_3+1$. In particular $H(0)=1$ since it is nonnegative and divides $-1$. By comparing degrees and constant terms, for integers $b,t,h$ we obtain
$$
 L_1=bH-(H-1)/X,\quad L_2=tH,\quad
 L_3=hH+(H-1)/X,\qquad b+t+h=m .
$$
Here $t\ge1$ and $b\ge1$ by nonnegativity. The leading coefficient of $L_3$ in degree $n$ gives $h\ge0$. Since $U$ is coprime to $PH(P)$ and divides the second actual part, $U\mid t$, hence $t\ge3$. The third actual part is $(hP+1)H(P)$, so $h=0$ is impossible. If $b=1$, nonnegativity of $L_1$ gives $1=H_0\ge H_1\ge\cdots\ge H_n\ge1$. Then $L_1=X^n$ and $U\mid P^{n+1}$, impossible. Thus $b\ge2$, $h\ge1$ and $t\ge3$. At $m=6$ equality throughout would give $U=3$, contradicting $\gcd(U,m)=1$. This proves condition 1. The same identities give $H_{i+1}\le bH_i$ and $b\le m-4$, proving conditions 2 and 3.

For condition 4, $U\le t<m$ so the other row cannot take the prime $r$ of $m=r^e$. Its three positive normalized parts of sum $U$ would need two multiples of $m$. Thus $A$ takes $r$. If $r\mid P$, then $H(P)$ is a unit at $r$, and the unshifted part $tPH(P)$ must be zero at the full $r$-power in $mPH(P)$, since the extra 1 is in slot 3. This would give $m\mid t$. If $r\mid H(P)$ then $r\nmid P$ and write the actual parts as
$$((bP-1)H(P)+1,\ tP H(P),\ (hP+1)H(P)).$$
Now slot 1 must be exceptional. Again the full power gives $m\mid t$. Hence $PH(P)$ must be a unit at $r$. The middle part is nonzero modulo $m$ because $0<t<m$, so it must be the exceptional slot. We deduce $hP\equiv-1\pmod m$ with $1\le h\le m-5$. This rules out $jP\equiv1\pmod m$ for $1\le j\le4$, as claimed.

For example, at $m=7$ the mixed case requires $h=1$ or 2, hence $P\equiv6$ or $3\pmod7$. Every square prime power and every power of two therefore gives the shifted conclusion for every $H$ satisfying the irreducibility hypotheses and the threshold.

## Unshifted weights

Here $U=Pv$, and $v$ divides $G(P)$ and all $L_i(P)$. Unless all three $L_i$ are multiples of $H$ as polynomials, $v\le R^{2n}$ by a nonzero resultant. As in the primitive case the assigned predecessor powers have product $S$ with $S^2\le R^{6n+3}$ by irreducibility of $XG-1$. Then
$$27(PG(P)-1)^2\le P^3 R^{12n+3},$$
contradicting $P>R^{15}$ since $PG(P)-1\ge P^{n+1}/2$ and $(12n+3)/(2n-1)\le15$.

In the remaining case $L_i=t_iH$ with positive integral $t_i$ summing to $m$. Pick a coefficient of $H$ coprime to $p$, which exists by primitivity. At this coefficient the bottom $v_p(m)$ base-$p$ digits of the total vanish. Carry-free addition requires them to vanish in every summand. Thus $p^{v_p(m)}$ divides every $t_i$, forcing $m/p^{v_p(m)}\ge3$. At any nonzero coefficient, three positive summands require at least three base-$p$ digit units. These prove both unshifted criteria.

Run `python3 check_polynomial_content.py` with SymPy to check exact polynomial division on 55 proportional and 57 mixed splits. The examples include 15 mixed specializations with row gcd at least three.
