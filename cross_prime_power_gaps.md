# Simultaneous leading prime powers

Let $A,C$ be any two positive rows of weight $N$, and set $D=\gcd(N,K_A,K_C)$ and $J=\gcd(N-1,K_A,K_C)$. Select two different positions in each row, with entries $a,a'$ and $b,b'$, respectively. Then
$$N^2(N-1)\mid D^2J\,aa'bb'. \tag{1}$$
In particular, in a counterexample every such pair of sums from the two rows other than the endpoint satisfies
$$(a+a')(b+b')\ge4N\sqrt{N-1}. \tag{2}$$
Indeed, in the [retained factor assignment](endpoint_pair_bound.md), $U^2V^2$ divides $aa'bb'$. Each power assigned from $N-1$ divides one of the two entries selected from its row, since only one entry can have residue 1. Thus $ST$ divides $aa'bb'$ and is coprime to $UV$. Use $DUV=N$, $JST=N-1$ for (1). The [endpoint theorem](all_weights_endpoint_reduction.md) gives $D=J=1$ in the counterexample. Two applications of AM–GM give (2).

The [two-row divisibility theorem](two_row_sums.md) strengthens the exponent in (2) to $20/13$ and gives intervals at every degree from five upwards.

## A criterion using two bases

Suppose $P=p^i$, $Q=q^j$ with $p,q$ distinct primes, $i,j\ge1$, and $PQ\mid N$. Choose powers $P_*=p^u$ and $Q_*=q^v$ in $(N/2,N]$. Put $h=N-P_*$, $k=N-Q_*$ and $g=\min(h,k)$. If
$$hk<4N\sqrt{N-1} \tag{3}$$
and for some $W=q^m$, $m\ge1$,
$$[-P]_W\,\left\lfloor\frac{g}{PQ}\right\rfloor
   <2W-\left[\frac NQ\right]_W, \tag{4}$$
the full conjecture holds at $N$. Brackets denote the least nonnegative residue. One may interchange the roles of $p,q$ in (4). A simple sufficient form of (4) is $[-P]_Q\lfloor g/(PQ)\rfloor\le Q$.

To prove the criterion, a row avoiding $p$ places the single leading $P_*$ in one part by carry freedom. Its two other parts have sum at most $h$. Similarly a row avoiding $q$ has two parts summing to at most $k$. By (2) and (3) the bases cannot be avoided by different rows of a counterexample paired with the endpoint. Since $p,q\mid N$, both must therefore be avoided by one row.

Every part of that row must be divisible by $PQ$. Its two smallest parts are $PQx,PQy$ with $x,y>0$ and $x+y\le\lfloor g/(PQ)\rfloor$. On dividing all parts by $Q$, their addition still has no base-$q$ carries. In particular $[Px]_W+[Py]_W\le[N/Q]_W$. Writing $s=[-P]_W>0$ gives $[Px]_W\ge W-sx$ and $[Py]_W\ge W-sy$, contradicting (4).

## An infinite family around fourth powers

Take any prime powers $P<Q$ of distinct bases, with $s=Q-P\le Q/10$. The full conjecture holds at **every multiple of $PQ$** in the interval
$$Q^4\ \le\ N\ \le\ Q^4+\frac{Q^3}{4s}. \tag{5}$$
There is no condition on the factorization of $N/(PQ)$ or $N-1$.

Indeed $P\ge0.9Q$ and $Q\ge10$, so $P^4,Q^4>N/2$. Write $k=N-Q^4$, $h=N-P^4$. Since $Q^4-P^4\le4sQ^3$,
$$hk\le Q^6(1+1/(16s^2))<4N\sqrt{N-1}.$$
Also $[-P]_Q=s$ and $sk\le Q^3/4<PQ^2$, proving (4).

For every $\epsilon>0$ and all sufficiently large $X$, **at least $\gg_\epsilon \sqrt X/\log X$ distinct weights $N\le X$** satisfy (5) with *prime* $P,Q$ and $N-Q^4\ge N^{3/4-\epsilon}$. Thus both leading powers can be at distances well above the $N^{4/7}$ scale for a single leading power.

To count them, take primes in $[Y,2Y]$. By the prime number theorem their number $m$ is $\gg Y/\log Y$. Let $F(v)$ count their unordered pairs at distance at most $v$. Partition the interval into $\lceil Y/v\rceil$ subintervals, with a shared endpoint assigned to only one. By Cauchy–Schwarz the pairs within the same subinterval number at least $(m^2/\lceil Y/v\rceil-m)/2$. Hence $F(v)\gg Yv/\log^2Y$ for $C\log Y\le v\le\sqrt Y$ and an absolute sufficiently large $C$.

Take $\theta=\min(2\epsilon,1/2)$ and $H=Y^\theta$. Summing $1/s$ over pairs with $s=Q-P\le H$ gives at least $\int_{C\log Y}^H F(v)\,dv/v^2\gg_\epsilon Y/\log Y$. For each pair use the upper half of (5). Its length in units of $PQ$ is $Q^2/(8Ps)\ge Y/(8s)$, so it contains $\gg Y/s$ multiples. Each has $N-Q^4\ge Q^3/(8s)\gg Y^{3-\theta}>N^{3/4-\epsilon}$ for all sufficiently large $Y$. These weights are at most $17Y^4$ for large $Y$, and each then has at most four distinct prime divisors in $[Y,2Y]$. Thus at most six counted pairs give the same $N$. Taking $Y=(X/17)^{1/4}$ proves the count.

## The other three residue classes

The same fourth powers give a certificate when the factors belong to either of the consecutive integers $N,N-1$. Assume $2\le s=Q-P\le Q/10$. The full conjecture holds at every weight satisfying
$$Q^4\le N\le Q^4+\frac{Q^3}{4s^2},\qquad
 N\bmod P\in\{0,1\},\qquad N\bmod Q\in\{0,1\}. \tag{6}$$
Thus all four residue classes modulo $PQ$ are covered in this shorter interval. The case $P,Q\mid N$ also has the longer interval (5).

As before (3) forces the same row to avoid both base primes. Its two smallest parts $a_1,a_2$ have sum at most $g=N-Q^4$. In a counterexample all parts of the nonendpoint rows are at least 3 since their scales are sums of three positive parts. For $a=a_i$, write $\delta=[a]_P$ and $\eta=[a]_Q$. Each is 0 or 1, and each of $\delta,\eta$ occurs with value 1 at most once in the row, by the low-digit test. The integer
$$z=\frac{sa-Q\delta+P\eta}{PQ}$$
is positive. The numerator is a multiple of $PQ$ greater than $-PQ$. If it were zero, the only possible positive values of $a$ would be $1$ and $Q/s$. The latter is not integral since $s>1$ and $\gcd(Q,s)=1$.

Let $d_P,d_Q$ be the least positive inverses of $s$ modulo $P,Q$. At either modulus $M=P,Q$ one has $\min(d_M,M-d_M)\ge(M-1)/s$, since both $s d_M$ and $s(M-d_M)$ are positive and congruent to $1$ and $-1$ modulo $M$, with $s>1$. The sum of the two values of $z$ satisfies
$$z_1+z_2\le\frac{sg+P}{PQ}\le\frac{Q^2}{4sP}+\frac1Q < \frac{P-1}{s}.$$
The last bound uses $P\ge0.9Q$ and $Q\ge20$.

Put $\sigma=\delta-\eta$. Direct division gives the next blocks of digits of each part:
$$\left[\left\lfloor a/Q\right\rfloor\right]_Q=[-z+d_Q\sigma]_Q,\qquad
  \left[\left\lfloor a/P\right\rfloor\right]_P=[z+d_P\sigma]_P.$$
If either $\sigma$ is zero, its $Q$-block is $Q-z$. The other $Q$-block is at least $\min(d_Q,Q-d_Q)-z'$ (and is $Q-z'$ when its sign is zero). Their sum exceeds $Q$ by the bound on $z_1+z_2$. If both signs are nonzero, one is positive and one negative, by uniqueness of the low digits. Their $P$-blocks are then $z_1+d_P$ and $z_2+P-d_P$ after reordering, with no wrap at either entry. Their sum exceeds $P$. Either case forces a carry. The proof in fact only needs (3) and $(sg+P)/(PQ)<\min(d_P,P-d_P,d_Q,Q-d_Q)$, with $g$ the smaller leading-power gap and $s\ge2$.

**Consecutive powers ($s=1$).** All four patterns in (6) work in this case as well. The same $z$ is now nonnegative. In the $Q$-block the two entries have residues $[-u_i]_Q$ with $u_i=z_i-\delta_i+\eta_i$. We have $u_1+u_2\le g/(PQ)+1/Q+1<Q$ at the stated bound. If both $u_i\ge1$ their residues already force a carry. If $u_i\le0$, the formula $a_i=PQz_i+Q\delta_i-P\eta_i$ leaves only $a_i=1,Q,Q^2$ among positive parts. The first is excluded by the row scale $U\ge3$. If $Q\mid N-1$, then $U$ is coprime to $Q$, also excluding the other two. The pattern $P,Q\mid N$ follows from (5).

It remains to treat $P\mid N-1$, $Q\mid N$. Here $N=Q^4+PQm$ with $0\le m\le Q^2/(4P)<Q$. For $m=0$ the weight itself is a prime power, so every positive row has a carry. For $m>0$ let $v$ be the full $q$-part of $m$. The full $q$-part of $N$ is $Qv$ since $v<Q$. As $U\mid a_i\in\{Q,Q^2\}$ and $U\ge3$, the row takes precisely this power from $N$. All the powers $S$ assigned to this row from $N-1$ must divide $a_i-1$, since $a_i$ is coprime to $N-1$. Hence
$$S\mid\gcd(N-1,Q^2-1)=P\gcd(m,Q+1),\qquad S\le Pm/v.$$
For the other row each power in $T=(N-1)/S$ divides two of its three normalized parts of sum $U$. Therefore $27(N-1)^2\le U^3 S^2\le Q^3 P^2v m^2<Q^8$, impossible for $N\ge Q^4$.

For each of the four fixed residue patterns in (6), at least $\gg \sqrt X/(\log X)^3$ distinct weights $N\le X$ satisfy that pattern and (6), with
$$N-Q^4\gg N^{3/4}/(\log N)^2.$$
Among the $\gg Y/\log Y$ consecutive prime pairs within $[Y,2Y]$, at least half have gap $s\le C\log Y$ for some absolute $C$, since the sum of the gaps is at most $Y$. These gaps are even for large $Y$. For each such pair the upper half of (6) contains $\gg Y/(\log Y)^2$ integers in the prescribed residue class modulo $PQ$. Each has the required distance. At most eight primes in $[Y,2Y]$ divide $N(N-1)$ for such a weight when $Y$ is large, so multiplicity is bounded. Again use $Y=(X/17)^{1/4}$.

## The full three-quarter scale near cubes

For prime powers $P<Q$ of distinct bases with $2\le s=Q-P\le Q^{1/4}/4$, the full conjecture holds at every
$$Q^3\le N\le Q^3+Q^{9/4},\qquad N\bmod P,\ N\bmod Q\ \in\{0,1\}. \tag{7}$$
Here both leading powers can be at distance at least $\tfrac12 N^{3/4}(1-o(1))$ from the weight. For each of the four fixed residue patterns there are $\gg\sqrt X/(\log X)^2$ such weights $N\le X$ in the upper halves of these intervals.

Indeed $Q\ge4096$ and $P>0.99Q$, so $P^3,Q^3\in(N/2,N]$. Put $g=N-Q^3$. Since $Q^3-P^3\le3sQ^2$, the product of the two gaps is at most $(7/4)Q^{9/2}<4N\sqrt{N-1}$. Also
$$s(sg+P)\le Q^{11/4}/16+sP <PQ(P-1).$$
This yields the mixed-block test proved above. For the quantitative assertion take primes $P,Q\in[Y,2Y]$ with $s\le Y^{1/4}/4$. The bound for $F(v)$ above gives $\gg Y^{5/4}/\log^2Y$ such pairs. For each, the upper half of (7) gives $\gg Y^{1/4}$ weights in each CRT class. All are at most $9Y^3$ for sufficiently large $Y$. The number of such pairs for one weight is bounded since $PQ\mid N(N-1)$ and every counted prime is at least $Y$. Taking $Y=(X/9)^{1/3}$ proves the count.

## Checks

Run `python3 check_cross_prime_power_gaps.py`. It checks the leading-digit and joint-row obstructions against exact binomial products through weight 510, including situations in which one row really avoids both bases. The small-weight run gives 242 joint-row obstructions and 111,717 separate-row comparisons. A separate next-block check tests 26,362 pairs of parts, including 3,856 exclusions by the mixed-block bound and 5,526 pairs with both next blocks compatible. The checker also verifies the consecutive-power calculation at 16,460 multipliers for four pairs of powers. The [cross-part verifier](check_cross_parts.py) checks (1) directly. Some examples of (5) are below. The endpoints in the table are the first and last certified multiples.

| $P$ | $Q$ | First $N$ | Last $N$ | Number of multiples |
|---:|---:|---:|---:|---:|
| 29 | 31 | 924,172 | 926,869 | 4 |
| 1,009 | 1,013 | 1,053,023,773,846 | 1,053,087,145,100 | 63 |
| $2^{19}$ | $3^{12}$ | 79,766,443,077,089,260,732,416 | 79,766,443,082,104,567,234,560 | 19 |
| 1,000,000,007 | 1,000,000,009 | `1000000036000000487000002924000006552` | `1000000036125000489000002931875006552` | 125,000,001 |
