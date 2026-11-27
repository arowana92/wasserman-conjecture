# Cofactor ten when one varying power divides the other minus one

Let $P<Q$ be powers of distinct primes, with $P\mid Q-1$. **The full trinomial conjecture holds at $N=10PQ+1$ whenever $Q$ has at most four nonzero digits in base $P$. This includes every $Q\le P^4$. If either power is even, the conclusion holds for every $Q$ without an upper bound.** A [separate residue argument](nested_power_residues.md) removes the digit condition whenever $P$ is a square or $Q$ is a fourth power. Both varying powers can be ordinary primes. For fixed even $P$ there are, in particular, infinitely many ordinary primes $Q$ in the indicated class.

## Determining the candidate row

Write $Q=mP+1$. By the [scale test at ten](static_square_cofactors.md), a possible failure requires $5\nmid PQ$ and $(U,x)=(3,3),(3,6),(3,9)$ or $(7,7)$. The nonendpoint row $A$ of scale $U$ takes 5 and the bases of $P,Q$, at three distinct exceptional positions as in the [digit-gap lemma](predecessor_inverse_digit_gaps.md). Its part exceptional at 5 is $xPQ$. Its part exceptional at $Q$ has the form $PL$ with $L\equiv-m\bmod Q$, since $mP=Q-1$. The full row is therefore
$$
 A=\big(xPQ,\ P(tQ-m),\ Q(kP+1)\big),\qquad x+t+k=10,\quad t,k\ge1. \tag{1}
$$
Indeed $t\ge1$ by positivity, $t\le10-x$ since $m<Q$, and $k=0$ would give a part $Q$ coprime to $U\ge3$. Both nonexceptional parts at 5 are multiples of $5U$. In particular
$$
 kP+1\equiv0\bmod 5U,\quad xPQ\equiv1\bmod5,\quad U\mid N,\qquad \gcd(U,N/U)=1. \tag{2}
$$
Here $U$ is prime, so the last condition excludes $U^2\mid N$. The scale $7$ permits a binary carry-free companion only in the shape $(V,2V,4V)$ with $V=N/7$: its normalized positive parts have sum seven and distinct 2-adic valuations, of minimum total $1+2+4$. At scale three the three identical companion parts force a binary carry, so $A$ must avoid 2.

## Either power is even

The full power of 2 in $N-1$ includes the extra factor 2 in 10. Hence the part $xPQ$ is divisible by twice the even varying power, so $x$ is even. This forces $x=6$, $U=3$, $k+t=4$, and $15\mid kP+1$. If $P$ is even, $k=3$ is impossible modulo 3; at $k=1,2$ the product $kP$ is a power of two, which cannot be $-1\bmod15$.

If $Q$ is even, $Q(kP+1)$ must vanish modulo $2Q$, so $k$ is odd. Thus $k=1$, $t=3$, and $P\equiv-1\bmod15$. Above the $Q$ shift the three binary parts
$$
 6P,\quad 3P-1,\quad P+1
$$
must have disjoint bits. In fact this requires $3P=2^{2j+3}+1$ for some $j\ge0$. Then $P\equiv-1\bmod5$ would give $2^{2j+3}\equiv1\bmod5$, impossible at an odd exponent.

For the binary assertion read the bits of odd $P$ from right to left, with initial multiplication carries $(0,-1,1)$ for the three affine expressions. On reading bit $e$, add $(6e,3e,e)$ to the carries. Reject when at least two resulting values are odd; otherwise replace the values by their floors after division by two. The first bit 1 gives state $(3,1,1)$. The complete reachable transition table is below, with a dash for rejection.

| State | Next bit 0 | Next bit 1 |
| --- | --- | --- |
| $(3,1,1)$ | – | $(4,2,1)$ |
| $(4,2,1)$ | $(2,1,0)$ | $(5,2,1)$ |
| $(5,2,1)$ | – | – |
| $(2,1,0)$ | $(1,0,0)$ | $(4,2,0)$ |
| $(4,2,0)$ | $(2,1,0)$ | – |
| $(1,0,0)$ | $(0,0,0)$ | – |
| $(0,0,0)$ | $(0,0,0)$ | – |

The low bits are $1,1,0$, followed by zero or more repetitions of $1,0$, then zeros. Summing their place values gives the assertion.

## Bounding the digits of the quotient

Suppose both powers are odd, $P>500$, and $Q$ has at most four nonzero base-$P$ digits. Write $m=\sum_{j=0}^{L-1}z_j P^j$, using ordinary base-$P$ digits. There are at most three nonzero $z_j$. With $z_{-1}=1$, we claim that $L\le3$ and
$$
 1\le z_j\le t z_{j-1}\quad(0\le j<L). \tag{3}
$$
Above the $P$ shift the row (1) adds without carries as
$$
 x+xPm,\quad t+(tP-1)m,\quad k+(kP+1)m .
$$
At position $j$ the raw digit coefficients are $xz_{j-1},t z_{j-1}-z_j,k z_{j-1}+z_j$. Assume the previous raw coefficients have been ordinary nonnegative digits, so none borrows or carries into position $j$. Consider the first three positions, up to the top digit if earlier. Since $k+t\le7$, we have $10k t^{j}<P$ for $j=0,1,2$: the maximum $k t^2$ is 50. In particular $10z_{j-1}<P$ by induction. If both $z_j>t z_{j-1}$ and $z_j<P-kz_{j-1}$, reducing the three coefficients modulo $P$ adds an extra $P$ to their sum, contradicting carry freedom.

The other possibility outside the claimed upper bound is $z_j=P-s$ with $1\le s\le k z_{j-1}$. Grouping the last two upper parts shows that $xQ,(10-x)Q$ add without carries. Their units blocks are $x,10-x$, so $xm,(10-x)m$ also add without carries. At digit $j$ these two parts have residues $P-xs,P-(10-x)s$, with no multiplication carries from earlier digits by $10z_i<P$ for $i<j$. Since $10s<P$, these residues already exceed $P$ in sum. Thus $0\le z_j\le t z_{j-1}$ and all three raw coefficients at $j$ lie between zero and $10z_{j-1}<P$, completing the induction. If the first three digits are nonzero then there are no higher nonzero digits. If a zero occurs among them, the previous argument leaves no multiplication borrow or carry above that position. At the next position the raw digits are $0,-z,z$, whose residues sum to $P$ whenever $z>0$. Thus zero propagates through all remaining positions without any size hypothesis on their exponents. This proves (3) and $L\le3$.

## Exhausting the polynomial shapes

Each digit tuple now fixes $m$ as a polynomial in $P$. Use (2) modulo $5U$ and the fact that $m$ is even at odd $P,Q$. For $L=1,2,3$ the numbers of tuples $(U,x,k,t,z_0,\ldots,z_{L-1})$ remaining are respectively **0, 9, 470**. All have $U=3$. Hence $A$ must avoid $2,3,5$. For a base $r$ one can retain just the residues of $P$ modulo $r^j$ passing the required congruences and having no carry through the first $j$ digits of (1). Lift each remaining residue by one digit. For all nine two-digit tuples and 469 of the three-digit tuples one of these sets is empty, at depths at most 60 in base two or 25 in base three. The verifier independently generates the shape list by enumeration of all static residue classes and rechecks every empty tree using the individual digits of the three polynomial parts.

The last tuple is
$$
 (U,x,k,t)=(3,6,2,2),\qquad m=2+4P+4P^2 .
$$
Put $T=2P+1$. The second part is $4P^3T$, and the sum of the second and third parts is $T^2(4P^2+1)$. Since $5\mid T$ and $5\nmid P$, the second part has valuation $v_5(T)$ whereas their sum has higher valuation. Thus the addition carries at 5. [Lean-checked integer identities](PredecessorPairWindows.lean)

## Exhausting the small powers at unrestricted digit positions

For $P\le500$ use each odd prime power $P$ other than powers of 5, each permitted $(U,x,k,t)$, and at most three occupied positions of $m$. The following stopping rule makes the search finite. At place $B=P^j$, suppose the lower digits give $r=m\bmod B$ with $10(1+Pr)<B$. At $m=r$ the three upper parts are nonnegative and each is smaller than $B$, by their sum. If there were a higher nonzero digit, take the first, of value $z$ at $B'\ge B$. At that position the upper parts for the full $m$ then have digits $0,P-z,z$, even allowing every more significant digit. This forces a carry. Consequently the prefix cannot be extended.

Start at units and prune whenever the sum of the upper-part residues crosses the current block boundary. Apply the stopping rule after each choice or stop after three nonzero digits. The numbers of positive quotients reaching a terminal test with one, two and three occupied positions are respectively 377, 3,794 and 142,245, counting the choice of row parameters. At each such quotient check (2) and oddness, then the actual factorial-floor valuations at 5, $U$, the base of $P$ and 2, with the binary companion possibility at $U=7$. No case survives.

An independent enumeration uses the [base-prime gap bound](predecessor_inverse_digit_gaps.md). For every permitted row parameter in this finite range, $P=p^a\ge10-x$. The gap between successive occupied base-$p$ positions in $Q$ is therefore at most $2a$. Successive occupied base-$P$ positions can differ by at most two, since a difference of three would force a base-$p$ gap of at least $2a+1$. This applies also above the units digit of $Q$, which gives the initial position bound for $m$. The second verifier places each nonzero digit by the available digit budget of the second part and tests the resulting rows by direct base-prime additions. It excludes all 3,400 candidates reaching its static congruence filter. Run **python3 check_ten_nested_powers.py** to reproduce the calculations. [Standard-library verifier](check_ten_nested_powers.py)
