# Modular inverses in the second interval above a power

Let $P=p^a$ and $Q=q^b$ be powers of distinct primes, with $2\le P\le Q$ and $PQ\mid N-1$. Put $m=PQ$, $r=(N-1)/m$, and let $L>N/2$ be a prime power at most $N$ whose base divides $N(N-1)$. The [two-power bound](two_predecessor_pair_sum.md) settles $N-L\le m$ and leaves at most one possible pair sum when $m<N-L\le2m$.

**Throughout $N-L\le2m$ the full conjecture holds if any of the following is true.**

* $P\le128$ and $P\ne101$, with no upper bound on $Q$ or $N$.
* $P\le128$ and $Q\le10^{1000}$. At $P=101$ the unrestricted result also holds unless $r=22$ and $Q$ is an odd power of a prime $2\bmod3$.
* $27r^2>(P-1)^3$.
* $r-1$ is prime. This permits arbitrarily many other factors in $r$.
* In the base-$p$ expansion of $Q$, two successive occupied digit positions differ by more than $a$. The symmetric test uses the base-$q$ digits of $P$ and the exponent $b$.

These tests also hold for the sum of **any** two parts of either nonendpoint row at a hypothetical counterexample: under any listed condition that sum exceeds $2m$. The statements for weights follow by taking the row avoiding the base of $L$ and the two parts missing $L$.

For **every $P\le128$**, the full conjecture holds throughout the second interval when $L$ is a power of two or a power of $q$, without an upper bound on $Q$. In particular for every such modulus $PQ$ the first two integers congruent to 1 modulo $PQ$ above each power of two exceeding $2PQ$ both satisfy the conjecture.

## Only one pair of parts

In the remaining case the two parts $x,y$ sum to $m+1$. Their row $A$ avoids both $p,q$, with one part a multiple of $P$ and the other a multiple of $Q$. Define
$$u=Q^{-1}\bmod P,\qquad v=P^{-1}\bmod Q$$
by the least positive residues. Up to order, necessarily
$$x=Pv=(P-u)Q+1,\qquad y=Qu,\qquad z=N-x-y=(r-1)PQ.$$
Indeed the positive quotients in $Pv+Qu=PQ+1$ are strictly below $Q,P$ respectively, and the two residues are fixed. In particular $\gamma=\gcd(x,y)=\gcd(u,v)$, coprime to $PQ$. The row scale $U$ is a unitary divisor of $N$ with $3\le U\mid\gamma$. Its primes must avoid carries even in $x+y$, since $\binom{x+y}{x}$ divides the trinomial coefficient. The same necessary condition holds at $p$ and $q$. One can check all of these with the two fixed parts alone. If $P=2$ or 3 then $u\le2$, which already rules out the pair.

## A digit gap independent of the weight

In base $p$ the pair adds freely only if $v$ and $Q-v$ add without a carry to $Q$, since $x=Pv$ and $y=P(Q-v)+1$. Put $s=Q-v=(uQ-1)/P$. Suppose an occupied digit of $Q$ at position $t$ has at least $a$ zeros immediately below it. Write $Q=H p^t+D$ with $p\nmid H$ and $0<D<B=p^{t-a}$. Positivity of $D$ follows from $p\nmid Q$. Then
$$s=u H B+(uD-1)/P,\qquad 0\le (uD-1)/P<B.$$
The digit of $s$ at position $t-a$ is $uH\bmod p\ne0$, whereas $Q$ has zero there. Equivalently $s\bmod(pB)>Q\bmod(pB)$, ruling out a carry-free sum. This last modular inequality is [checked in Lean](InverseDigitGap.lean) under the displayed integer hypotheses.

For fixed $P>2$, the number of integers $Q\le X$ with $p\nmid Q$ that pass the necessary gap test is $O_P(X^\theta)$ where $\theta=\log_P(P-1)<1$. Indeed all digits in the base-$P$ expansion must be nonzero to pass the test. An internal zero block is followed (possibly after other zeros) by a nonzero base-$p$ digit. There are at most $(P-1)^k$ choices at length $k$. In particular all but $O_P(X^\theta)$ primes $Q\le X$ satisfy the second-interval conclusion simultaneously for all eligible $N,L$.

For fixed distinct primes $p,q$ and fixed $a$ with $P=p^a>2$, among exponents $1\le b\le B_0$ at most $O_{P,q}(B_0^\theta)$ can fail to give the same conclusion. Choose $k=\lfloor\log_P B_0\rfloor$. Apart from $O_{P,q}(\log B_0)$ exponents with $q^b<P^k$, the low $k$ base-$P$ digits must all be nonzero. The [order of $q$ modulo $P^k$](prime_power_ray_density.md#counting-exponents) is $\gg_{p,q}P^k$, so each residue accounts for $O_{p,q,P}(1)$ exponents. This proves the estimate; the case $P=2$ is already fully excluded by $u\le2$.

## Bounding the cofactor without the other power

Assign the full powers of $N$ as $UV=N$, and those of $N-1$ as $ST=N-1$, with $U,S$ assigned to $A$. Since $A$ takes $p,q$ at opposite exceptional positions, their full powers divide $x+y-1=m$. Thus $P,Q$ are already full in $N-1$ and $\gcd(r,PQ)=1$. Bins with exceptional position $x$ or $y$ contribute exactly $m$: their product divides $x+y-1$ and contains $PQ$. The remaining bin $w$ divides both $x/U$ and $y/U$. Therefore
$$S=mw,\quad r=wT,\quad Uw\mid\gamma,\quad r\equiv1\bmod U.$$
The last congruence uses $U\mid N$ and $m\equiv-1\bmod U$. The companion row equals $V(f,g,h)$, with $f+g+h=U$, and $T^2\mid fgh$. AM–GM gives
$$27r^2\le w^2 U^3\le\gamma^2 U\le (P-1)^3.$$
Here $r>1$ because $z$ is positive. For each fixed smaller power only finitely many cofactors, scales and normalized companion triples can occur, even as $Q$ runs through unbounded prime powers.

There is another consequence of $U\mid r-1$. If $r-1=s$ is prime, then $U=s$ (or $U\ge3$ is already impossible). Write $h=x+y=PQ+1$, and $s^e\parallel h$. The final part is $z=s(h-1)\equiv-s\bmod s^{e+1}$, while the residue of $h$ at this modulus is at least $s^e\ge s$. Adding $h$ and $z$ therefore carries at $s$, which $A$ must avoid. More generally, if $r-1=s^j$ then the same argument requires $1\le v_s(h)<j$. It also excludes $r=5$, as the remaining valuation one cannot support $U\ge3$.

In particular the [exponent-window result](two_predecessor_pair_sum.md#entire-rays-of-exponents) with $r-2$ in the denominator applies at every $r\ge4$ one above a prime, for arbitrary distinct $p,q$. For example it applies to $1+462\cdot5^a13^b$, with six distinct primes in the predecessor. Since $461$ is prime, every $a$ works whenever $460\cdot13^b<5^k\le462\cdot13^b$ for an integer $k$.

## The classification through 128

For each $P$ enumerate $3\le u<P$ coprime to $P$, and divisors $U\ge3$ of $u$. Enumerate $r>1$ congruent to 1 modulo $U$, with $27r^2\le u^2U$ and $\gcd(r,P)=1$. Of the full powers dividing $r$, any that do not divide $u/U$ are forced into $T$. Each must divide two of $f,g,h$. A reverse enumeration starts with $f+g+h=U$, takes $T$ dividing the lcm of the three pairwise gcds, and $w\mid u/U$ coprime to $T$. Set $r=wT$ and test $r>1$, $r\equiv1\bmod U$. These two procedures give identical candidate triples at every $(u,U)$.

At a candidate the rows as functions of the unknown $Q$ are given exactly by
$$A=(uQ,(r-1)PQ,(P-u)Q+1),\qquad C=(f,g,h)(rP Q+1)/U.$$
Apply the [shifted-row carry tests](shifted_directions.md) to this specified pair. They impose the low residues at powers of $rP$, the weight scale $U$, binary avoidance by at least one row, and existence of a carry-free input at each assigned base. For example, for slopes $a_i$ and offsets $b_i$ the digit transition on input $d$ takes each carry $c_i$ to $\lfloor(c_i+a_i d)/p\rfloor$ and rejects if the three residues sum to at least $p$. A breadth-first search from the prescribed offsets decides existence at that base; the states are bounded because the slopes are fixed. The tests allow different inputs at different primes, so rejection suffices.

Of 1,983 companion choices, only the following possibilities survive. Both require odd $Q$.

| $P$ | $u$ | $U$ | $r$ | $(f,g,h)$ | $Q$ congruence | Bases $A$ must avoid |
|---:|---:|---:|---:|:---|:---|:---|
|101|77|7|22|$(1,2,4)$|$8101\bmod15554$|$7,11,101$|
|127|104|4|13|$(1,1,2)$|$2297\bmod6604$|$2,13,127$|

Through $Q=10^{1000}$, simultaneous carry tests at the three required bases return no inputs on either progression. One implementation intersects digit-graph prefixes, and a separate search uses intervals on which a factorial-floor carry is forced. The [complete cofactor-1651 argument](cofactor_1651.md) excludes the second shape for every prime power $Q$ and excludes squares and $Q\equiv1\bmod3$ at cofactor 2222.

For the remaining $P=101$ shape with $L$ binary, $A$ must avoid 2, whereas $Pv$ and $Qu$ are both odd. If $L=q^k$, the offset bounds give $2020<q^{k-b}\le2121$. Enumerate the prime powers in this interval and the successive powers of each prime modulo 15554. Only $q=2081$ has $8101$ in its orbit, and it requires $b=214\bmod300$. Such a $Q$ is a square. This proves the unrestricted result for both markers.

## Reproduction

Run `python3 check_inverse_digit_gap.py` for small exact binomial controls and the digit-gap test. With python-flint installed, run `python3 check_inverse_pair_scales.py` for the two enumerations and both independent 1,000-digit checks. The enumeration of normalized companions is independent of the size of the larger prime power. Run `lean InverseDigitGap.lean` with Lean 4.12.0 to check the modular digit lemma.
