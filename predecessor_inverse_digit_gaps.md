# Gaps in the digits of a second predecessor power

Let $P=p^a$, $Q=q^b$ be powers of distinct primes, with positive exponents. Suppose a trinomial row at $N=DPQ+1$ avoids both bases and the positions congruent to one modulo $P,Q$ are different. Denote the part at the third position by $xPQ$, and put $d=D-x$. Then $x\ge1$, $d\ge1$. For every layer $H=p^j$ with $j\ge1$, write $Q=mH+h$, $0<h<H$. When $m\ge1$ and $dh<H$, the integers below exist:
$$
 P I-Q J=m,\qquad |I|\le dmh,\qquad |J|\le m(dPh+1)/Q. \tag{1}
$$
In particular $\|mP^{-1}\|_Q\le dmh$, where the norm denotes distance to the nearest multiple of $Q$ after choosing an inverse. If $p\nmid m$, then $H\le dPh+1$. Two successive occupied base-$p$ digit positions of $Q$ differ by at most
$$
 a+\lceil\log_p d\rceil. \tag{2}
$$
The same result holds on exchanging the powers. Both inequalities concern the powers actually present in the predecessor and allow arbitrary exponents.

Indeed the part exceptional modulo $Q$ is $PL$, with $PL\equiv1\bmod Q$ and $L>0$. The third summand after removing $xPQ$ is $P(dQ-L)+1$. Above the $P$ shift, $L$ and $dQ-L$ add without carries at $p$. Therefore $L=lH+c$ with $0\le l\le dm$, $0\le c\le dh$. Put $I=mc-lh$; since $mL=lQ+I$ we have $P I\equiv m\bmod Q$. This gives (1). If $p\nmid m$ then $J\ne0$, whence $Q\le m(dPh+1)$. Together with $Q=mH+h$ this gives $H<dPh+1$. If $dh\ge H$ that last inequality already holds. At an occupied position $v$ with previous occupied position $u$, choose $H=p^v$. Here $p\nmid m$ and $h<p^{u+1}$. Hence $p^{v-u}<p^{a+1}d$, giving (2).

For $H=P$, multiplication by $P$ in $PL\equiv1\bmod Q$ and substitution of $mP\equiv-h\bmod Q$ give
$$
 Q\mid l h^2-m c h-m^2 .
$$
If $Q>dm h^2+m^2$ and $dh<P$, the expression vanishes, forcing $h\mid m^2$. In particular for *ordinary prime* $Q$ this forces $h=1$, since $m<Q$ and $(m,Q)=1$. For a general prime power $Q$ every prime factor of $h$ must divide $Q$. This gives a test with no hypothesis on the prime base or exponent of $P$.

## Consequences for the full conjecture

For $D=b_0$ a prime power or $D=2b_0$ with $b_0$ an odd prime power, every possible failure must satisfy the digit-gap tests, with $d\le D-3$. By the [scale and numerator theorem](static_square_cofactors.md), the two varying bases go to one nonendpoint row of scale $U\ge3$, the static base of $b_0$ also goes there, and neither varying base divides $b_0$. The three exceptional positions at these bases are distinct: otherwise some part would be divisible by $b_0PQ$ and by $U$, too large for the weight. The part exceptional at $b_0$ is $xPQ$ and $U\mid x$. This proves the assertion even when a varying power is even. At cofactor ten, the alternatives $x=3,6,9,7$ give $d=7,4,1,3$ respectively. In the $x=9$ case the gap bound is just $a$, at both powers with the corresponding exponent.

As an example, for ordinary primes $P<Q=mP+h$, $1<h<P$, the conditions $7h<P$ and $7m h^2+m^2<Q$ certify $1+10PQ$ at every weight. The complementary residue $h=1$ has a [stronger separate result](ten_nested_powers.md).

Run **python3 check_predecessor_pair_gaps.py** for direct factorial-valuation checks of (1), the gap bound and the quadratic congruence. [Verifier](check_predecessor_pair_gaps.py)
