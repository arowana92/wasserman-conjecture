# A power saving for every positive quadratic at a prime power

Fix positive integers $M,h$. Among prime powers $P\le X$, the number for which Wasserman's full three-row conjecture can fail at
$$
                         N=MP^2+hP+1
$$
is $O_{M,h}(X^{(5+\beta)/6})$, where $\beta=\log_2((1+\sqrt5)/2)$ and $(5+\beta)/6<0.950$. Thus the conjecture holds for a proportion tending to one of the prime inputs, by the prime number theorem. The other factor $MP+h$ of the predecessor can have arbitrarily many prime divisors.

More precisely, put $\gamma=(1+\beta)/2$. In every interval of $K$ consecutive input integers in $[X,2X]$, the bound is
$$
              O_{M,h}\big(X^{2/3}+K^\gamma X^{2(1-\gamma)/3}\big).       \tag{1}
$$
In particular the number is $o(K)$ uniformly when $K/X^{2/3}\longrightarrow\infty$. If $h^2-4M$ is negative or nonsquare, the bound improves to $O_{M,h}(K^\gamma)$ in every such interval. Earlier pointwise results cover several combinations of coefficients and discriminants. The counting argument here includes reducible quadratics. [Quadratic specialization](quadratic_powers.md) · [Counting with an additional large prime power in the predecessor](polynomial_predecessor_pairs.md)

## The first row

Write $P=p^k$, $H=MP+h$ and $N=1+PH$. We can omit bounded $P$, so assume $P>\max(M,h)$. A counterexample contains the endpoint row. Distribute the primes of $N(N-1)$ between the two other rows, assigning each to a row avoiding it. Name the row assigned $p$ the first, and write $U,V$ for the products of full weight powers assigned to the first and second rows. Then $UV=N$, $\gcd(U,V)=1$, both scales are at least 3, and each divides all parts of its row. One of these rows must have odd coefficient. [Endpoint reduction](all_weights_endpoint_reduction.md) · [Prime assignments](endpoint_reduction.md)

At base $p$, the first row splits the two blocks of $H$. Its two unshifted parts are $P(aP+b),P(cP+e)$, for nonnegative coefficients bounded solely by $M,h$, and its third part is determined by the total. The two unshifted forms are nonzero. The binary count in the next section shows that any of these finitely many first rows can be odd at just $O_{M,h}(K^\gamma)$ inputs in an interval of length $K$. It is enough henceforth to count an odd companion.

If $U$ is bounded in terms of $M,h$, the companion has parts $(N/U)(f,g,l)$ with $f+g+l=U$ positive. There are finitely many such forms. Restrict to each class modulo $U$ on which $N/U$ is integral. The binary count for its first two parts again gives $O_{M,h}(K^\gamma)$.

Otherwise both unshifted linear forms must share the same polynomial factor of $N$. Indeed $U\mid N,aP+b,cP+e$. A nonzero constant form bounds $U$, as does $ae-bc$ when nonzero. For $a>0$, reducing $a^2N$ modulo $U$ shows $U\mid a^2-hab+Mb^2$. If this vanishes, $b>0$ and $a/b$ is an integer root $u$ of $Y^2-hY+M$. Set $t=h-u>0$, so $M=tu$, $h=t+u$. The exact row is
$$
   N=TW,\quad T=1+tP,\quad W=1+uP,\qquad
   W(Pi,Pj,1+(t-i-j)P),                                    \tag{2}
$$
where $i=b$, $j=e$ are positive and $i+j\le t$. This proves the nonsquare assertion once the binary count is established.

## Counting a quadratic cut in binary

We use the following elementary estimate. Take two integral polynomials of degree at most two in an integer index $v$, and their sum. Suppose all three values are nonnegative and their derivatives are positive on the interval considered. For each of the three affine derivatives let $e_i$ be the minimum of the 2-adic valuations of its slope and intercept, with $v_2(0)=+\infty$. Write $e=\max_i e_i$. The number of indices in any interval of length $K$ where the two values add without a binary carry is
$$
                           O\big(\varphi^e K^\gamma\big),
                    \qquad \varphi=(1+\sqrt5)/2 .             \tag{3}
$$
The implied constant is absolute. The same estimate applies when a common power of two divides all nonconstant coefficients: first divide the two polynomials by that power and round their constants down separately. Carry freedom of the initial pair implies carry freedom of the new pair. Use the derivatives of the new pair and its sum.

Translate the index to start at zero; this preserves the minima defining $e_i$. Put $m=\lfloor(\log_2 K)/2\rfloor$, treating $K<4$ directly. Write $v=v_0+2^m y$, $0\le v_0<2^m$. In the next $m$ digits of each part the quadratic Taylor remainder disappears modulo $2^{2m}$. The digits at these positions are therefore the low $m$ digits of an affine function of $y$ with slope equal to the derivative at $v_0$. The fixed offsets are nonnegative. If the lower digits already carry, discard this prefix.

Let $b$ and $w$ be the minimum and maximum valuations of the two slopes and their sum at this prefix. If $w<m$, divide out the common delay $b$, rounding offsets down. The remaining $m-b$ columns have delay spread $w-b$. Use the finite-window version of the [binary two-summand estimate](binary_affine_counts.md). For each of the $m-b$ input steps, retain the determined future delayed digit even if its position is beyond the tested window. For two summands there is just one delayed stream. Such a retained digit cannot collide with another future delayed digit, as the stream has only one value per column. In every current column the two undelayed digits and the recorded digit must sum to one. With $j$ pending ones give the prefix weight $\varphi^{-j}$. If the current digit is one, at most one choice passes and its relative weight is at most $\varphi$. If it is zero, the new delayed digits differ, and the relative weight sum is at most $1+1/\varphi=\varphi$. Terminal weights are at least $\varphi^{-(w-b)}$. Thus at most $\varphi^{m-b+w-b}$ choices of the tested input bits pass. The $b$ remaining input bits cost $2^b$, so the total is at most $\varphi^{m+w}$ since $2\le\varphi^2$. For $w\ge m$ the trivial $2^m\le\varphi^{2m}$ suffices. Each input residue modulo $2^{2m}$ occurs at most four times.

For an affine derivative on the $2^m$ prefixes, removing its content $2^{e_i}$ leaves either a constant odd residue modulo 2 at every index or an odd slope. In the latter case each class modulo $2^k$, $k\le m$, is taken equally often. Since $\varphi<2$, the sum of $\varphi^{\min(m,v_2)}$ over the prefixes is $O(2^m\varphi^{e_i})$. Summing over the three derivatives bounds the maximum in the same way. This gives (3) because $(2\varphi)^m\le K^\gamma$.

## The smallest group of predecessor powers

Return to (2). Set $B=\gcd(T,U)$; then $T/B\mid V$. Since $U\mid W\gcd(i,j)$, every prime of $U$ contributes its full exponent in $TW$ to $U$, and hence $B\mid\gcd(i,j)$. Put $S$ equal to the product of full powers in $H$ assigned to the first row. Each power with base different from $p$ divides $i$ or $j$: its two zero positions intersect the two unshifted parts and both $P,W$ are units there. For the base $p$ use the full predecessor power, then divide by $P$; the $p$-part of $H$ divides both $i,j$. Thus $S\mid\operatorname{lcm}(i,j)$. Both $B$ and $S$ have only boundedly many choices.

All full powers of $H/S$ are assigned to the companion. Partition them into three pairwise coprime groups by the position of the unique residue one. The smallest group product $d$ satisfies $1\le d\le(H/S)^{1/3}$. The part at that position is
$$
                Y=\frac{T}{B}\,\frac{H}{Sd}\, n,\qquad n\ge1 .
$$
Indeed it is divisible by the two relatively prime factors displayed. Since $Y<N$, we have $Hn<BSd W$. Using $tW=H-u$ gives $tn<BSd$. Put $L=BSd$, $E=L-tn>0$, so $n=O_{M,h}(d)$. The residue one modulo $d$ gives
$$
       (H/(Sd))n\equiv BW\pmod d,\qquad
       t n\,H/(Sd)\equiv-Bu\pmod d .                         \tag{4}
$$
Here $T^{-1}\equiv W\bmod d$, as $TW=N\equiv1$. In particular $\gcd(n,d)=1$. These statements include $d=1$.

For fixed $B,S,d,n,t,u$, (4) gives at most $t$ possible residues for $H/(Sd)$ modulo $d$. Thus $H$ lies in at most $t$ classes modulo $Sd^2$. Since $H=MP+h$, pulling back each class costs at most $M$ classes of $P$ modulo $Sd^2$. Refine to classes modulo $Q=BSd^2$, at a further bounded cost. Retain classes on which $T/B$ is integral. The two polynomials $Y,N-Y$ in the index $v$ of $P=P_0+Qv$ have integral coefficients. The companion being odd implies carry freedom in this cut.

## The derivative bound on a class

Put $m_d=v_2(d)$, $m_n=v_2(n)$. Divide the cut by $2^{m_d}$ and round constants down as allowed in (3). Before this division its three derivatives in the index $v$ are $d$ times
$$
       n(2tH-t^2),\qquad E(2H-h)-Mn,\qquad L(2H-h).           \tag{5}
$$
The middle expression is the third minus the first. Since $E\ge1$, $n=O(d)$, $d=O(X^{1/3})$ and $P\ge X$, both parts and both derivatives are positive for large $X$. The positive-degree coefficients are divisible by $2^{m_d}$, as also follows by expanding the term of degree two using $Q/L=d$. After division the extra factor $d/2^{m_d}$ is odd.

The largest derivative content exponent in (3) is at most $\max(m_d,m_n)+O_{M,h}(1)$. To verify this, first suppose $m_d$ exceeds a sufficiently large constant depending on $M,h$. Both $H$ and $L$ are divisible by $2^{m_d}$ and $n$ is odd. The first expression of (5) then has valuation exactly $v_2(t^2)$, the last exactly $v_2(L)+v_2(h)$, and the middle the smaller of these. This holds throughout the class. If $m_d$ is bounded, then so are $v_2(Q)$ and $v_2(L)$. The slope in $v$ of the first expression has valuation $m_n+O_{M,h}(1)$ and the last has bounded slope valuation. For the middle expression, either $v_2(E)$ is bounded, in which case so is its slope valuation, or $v_2(E)$ exceeds a sufficiently large constant. In the latter case $E+tn=L$ forces $v_2(n)$ bounded, and $E(2H-h)-Mn$ has the bounded valuation of $Mn$. This proves the claim.

Consequently the count for the class is
$$
                O_{M,h}\big((\varphi^{m_d}+\varphi^{m_n})
                              (1+K/d^2)^\gamma\big).
$$
The interval may contain a single index. The estimate covers that case too. For all $Z\ge1$, $\sum_{j\le Z}\varphi^{v_2(j)}=O(Z)$ by summing a geometric series. Sum over $n=O(d)$ and then over dyadic blocks $d\asymp J$ up to $O(X^{1/3})$. Each block costs $O_{M,h}(J^2(1+K/J^2)^\gamma)$. Both resulting powers of $J$ have positive exponent, giving (1). Taking $K$ comparable to $X$ and summing the dyadic input intervals proves the global theorem.

The [arithmetic checker](check_quadratic_prime_power_density.py) tests (4) against actual allocations of full powers, checks (5) by finite differences and compares the finite-window binary counts with direct digit addition. It examines 1,043,259 scaled rows and 109,120 upper binary words. Lean checks the congruence in (4), both quadratic differences behind (5), and the Taylor congruence used in (3). [Formal algebraic identities](QuadraticCompanionCut.lean)
