# A digit test for two prime powers at any weight

Let $W<Z$ be powers of different primes and let $WZ\mid n$. Write $u=Z-W$ and $n/W=\sum_{j=0}^e D_j W^j$ in ordinary base $W$ with $D_e>0$. Here $e\ge1$. Put
$$L=e\left(1+u+\sum_{j=1}^e D_j(1+u)^j\right).$$
Suppose $2L+2<Z$, and that $n/Z=\sum_{j=0}^e t_j Z^j$ with $2|t_j|+2<Z$ for all $j$ and at most one $t_j\le0$. **No positive binomial coefficient of weight $n$ avoids both bases.** The primes are consequently incompatible for positive trinomials at both $n$ and $n+1$. [Transfer proof](clustered_predecessors.md)

## Removing the unrestricted lowest digit

In a simultaneously avoiding binomial row, each part is $WC(W)$ with $C$ a digit polynomial of degree at most $e$, each coefficient bounded by the corresponding $D_j$. Since $WZ$ divides the part, $C(-u)$ is a multiple of $Z$. Its nonconstant terms have total absolute value at most $L$, and $0\le C(0)<W<Z$, so $C(-u)\in\{0,Z\}$. Write $\delta=C(-u)/Z$ and put $A(X)=C(X)-C(0)+\delta X$. Then
$$C(X)+\delta(X-W)=A(X)-A(-u)=(X+u)Q(X),$$
$$Q(X)=\sum_{j=1}^e A_j\sum_{i=0}^{j-1}X^{j-1-i}(-u)^i.$$
Thus the part divided by $Z$ equals $G(Z)$ for $G(Y)=(Y-u)Q(Y-u)$. For an integer polynomial let its coefficient norm be the sum of the absolute values of its coefficients. In each inner summand the norm of $(Y-u)^{j-i}(-u)^i$ is $(1+u)^{j-i}u^i\le(1+u)^j$. Therefore the norm of $G$ is at most $e\sum_{j\ge1} A_j(1+u)^j\le L$.

As $G(u)=0$ and $G(Z)>0$, each part polynomial has a strictly negative coefficient. The two parts cannot have negative coefficients in the same position and neither can have one opposite a positive $t_j$. Indeed, for a signed coefficient $c$ with $2|c|+2<Z$, the incoming borrow is $-1$ or $0$. If $c<0$ its normalized digit is greater than $Z/2$, and if $c>0$ it is less than $Z/2$. Addition without carries makes the two required negative positions distinct and leaves only one position available. This proves the test. [Lean formalization of borrow propagation and the two-position contradiction](SignedDigitCarries.lean)

## An upper bound for the number of positive parts

There is a finer test that uses only the ordinary digits of $n/W$ and $n/Z$. It applies to rows of any length. Allow the correction $\delta=C(-u)/Z$ to be an arbitrary integer this time. It lies between
$$l=-\left\lfloor\frac{\sum_{j\ {\rm odd}}D_j u^j}{Z}\right\rfloor,\qquad r=\left\lfloor\frac{D_0+\sum_{j\ge2,\ j\ {\rm even}}D_j u^j}{Z}\right\rfloor.$$
In the argument above the new polynomial $A$ has $A_1\in[l,D_1+r]$ and $0\le A_j\le D_j$ for $j\ge2$. In every case its associated part polynomial $G$ still has the positive root $u$ and takes a positive value at $Z$.

For $0\le k\le e$ the coefficient of $Y^k$ in $G$ is the following linear form in those bounded integers:
$$G_k=\sum_{j=\max(1,k)}^e A_j(-u)^{j-k}
\begin{cases}
j,& k=0,\\
\binom{j+1}{k+1},& k>0 .
\end{cases}$$
Indeed, the contribution from the monomial of degree $j$ in $A$ is $(Y-u)((Y-u)^j-(-u)^j)/Y$. Obtain lower and upper bounds $l_k,r_k$ for the form by taking the lesser and greater endpoint contribution in each summand. Suppose $1-Z<l_k$ and $r_k<Z$ at every position, and write $E_k$ for the ordinary digits of $n/Z$ in base $Z$, padded by zeros. Then the length of any positive row avoiding both bases is at most
$$\sum_{\substack{0\le k\le e\\l_k<0}}\left\lfloor\frac{E_k}{Z+l_k-1}\right\rfloor. \tag{2}$$
The incoming borrow stays in $\{-1,0\}$ under these bounds. Every strictly negative coefficient costs at least $Z+l_k-1$ of the available digit $E_k$. Each part needs at least one negative coefficient, giving (2). In particular a bound of one transfers from binomials at $n$ to trinomials at both weights. A bound of two directly forbids positive trinomial rows at $n$.

## Any fixed number of nearby prime powers

Let $m\ge3$ and $0=h_0<h_1=a<h_2=b<\cdots<h_{m-1}=h$, where $h=b$ for $m=3$. Suppose the $P+h_i$ are powers of distinct primes and set $F=\prod_i(P+h_i)$, $s=b-a$ and $e=m-1$. Assume
$$a\sum_{i=2}^{m-1}\frac1{h_i-a}<1,\qquad
(b+s)\sum_{i=3}^{m-1}\frac1{h_i-b}<1. \tag{1}$$
The empty sum is zero. For every integer $d\ge1$ with
$P\ge10e\,3^e d B,\qquad B=b^2\prod_{i=3}^{m-1}h_i,$
the bases of $P,P+a,P+b$ certify any two positive trinomial rows at **both** $dF$ and $dF+1$. The endpoint theorem proves the conjecture at both weights. For four factors the [quartic test](quartic_prime_certificates.md) also handles additional gap arrangements and gives a smaller constant.

To apply the digit test to each pair of these three powers, write $D_W(X)=d\prod_{P+h_i\ne W}(X+P+h_i-W)$, so that $dF/W=D_W(W)$. At $W=P$ the coefficients of $D_W$ are positive. At $W=P+a$ write $D_W=(X-a)H(X)$ with $H=d\prod_{i\ge2}(X+h_i-a)$. Its constant is negative and its other coefficients are strictly positive by (1). In fact for a product $H$ with positive shifts $g_i$ one has $H_j\le H_{j-1}\sum_i 1/g_i$ for $j\ge1$, by the elementary symmetric sums of the reciprocals. A single borrow normalizes the constant, so all higher digits are bounded by the absolute values of the corresponding polynomial coefficients. At any of the three powers $W$, the two differences from the other two of the first three powers have magnitude at most $b$. Each remaining difference is at most $h_i$. Consequently the sum of the absolute coefficients weighted by $(1+u)^j$ is at most
$$d(1+u+b)^2\prod_{i\ge3}(1+u+h_i)\le d\,3^e B.$$
In particular all the coefficients have magnitude far below $P$ and the digit descriptions at the smaller powers apply. The bound on the test parameter is $L\le e(1+b+d\,3^e B)\le2ed\,3^e B$.

At $Z=P+a$ only the constant of $D_Z$ is nonpositive. At $Z=P+b$ write $D_Z=(X-b)(X-s)H(X)$ with $H=d\prod_{i\ge3}(X+h_i-b)$. Its constant is positive. In degrees $2\le j\le e$, the coefficient $H_{j-2}-(b+s)H_{j-1}+bsH_j$ is strictly positive by (1), with missing coefficients of $H$ interpreted as zero. Thus only the linear coefficient can be nonpositive. Every coefficient is bounded in absolute value by $d\,3^eB$, so the large-base conditions in the test hold.

## Counting products with unequal gaps

For each fixed $m\ge3$ and $0<\epsilon<1/m$, there are
$$\gg_{m,\epsilon}\frac{X^{2/m}}{D(\log X)^3}$$
distinct products $F\le X$ where the same three primes work at **all** multipliers $1\le d\le D$, uniformly in integers $D\le X^{1/m-\epsilon}$. Each of the $m-3$ outer gaps contributes a logarithmic factor to the count because their product, together with $b^2$, controls the carry bound.

Put $Y=X^{1/m}/2$ and $T=Y/(10e\,3^eD)$. For $m\ge4$ choose scales $H_3<\cdots<H_{m-1}$, each a power of two, with $\log_T H_j$ in a separate fixed closed interval strictly above $1/e$. Take these intervals in increasing order, sufficiently close to $1/e$ that their upper endpoints sum to less than one. Put $H_2=(T/\prod_{j\ge3}H_j)^{1/2}$. All scales grow as positive powers of $Y$, uniformly in $D$. Successive scale ratios grow as positive powers of $Y$, including $H_3/H_2$ since $\log_T H_2<1/e$. For $m=3$ use the sole scale $H_2=\sqrt T$.

At a fixed choice of scales select primes from intervals nested from size $H_{m-1}$ down to $H_2$. In an interval of length at most $H_j$ call a set of $k$ primes sufficiently large at that level if $k\ge c_j H_j/\log Y$. The constants $c_j>0$ below depend only on $m,\epsilon$. Partition $[Y,2Y]$ into $O(Y/H_{m-1})$ intervals of the largest size. By the prime number theorem, for small $c_{m-1}$ the sufficiently large sets in them contain $\gg Y/\log Y$ primes. At each level, every upper-quarter prime in such a set is separated from all lower-half primes by more than $\delta_jH_j$ for a fixed $\delta_j>0$. In fact an interval of length $\delta_j H_j$ contains fewer than $k/8$ primes if $\delta_j$ is small, since $\log H_j\gg_{m,\epsilon}\log Y$. [Montgomery and Vaughan, equation (1.12)](https://personal.science.psu.edu/rcv4/personal/Publications/large_sieve.pdf)

At level $j\ge3$ choose the prime $P+h_j$ from the upper quarter. Distribute the lower half into $O(H_j/H_{j-1})$ cells of size at most $H_{j-1}$. For sufficiently small $c_{j-1}$, sufficiently large child sets contain $\gg k$ of these primes, and the construction can recurse with each such set. At level 2 choose $P+b$ from the upper quarter. Among lower-half primes take $P,P+a$ in the same cell of length at most $\delta_2 H_2/(2m)$. There are only $O_{m,\epsilon}(1)$ such cells, so one supplies $\gg k^2$ choices, all with $s>(m-2)a$.

Inductively a sufficiently large level-$j$ set supplies
$$\gg_{m,\epsilon} k\,\frac{H_2^2\prod_{i=3}^j H_i}{(\log Y)^j}$$
choices of $j+1$ primes. The bottom level gives $\gg k^3$. At each higher level the child sizes sum to $\gg k$, and the independent $\gg k$ upper-quarter choices supply the additional $H_j/\log Y$. Thus a full choice of scales gives $\gg_{m,\epsilon}YT/(\log Y)^m$ prime sets. They have $b\le H_2$ and $\delta_j H_j<h_j\le H_j$ for $j\ge3$. Hence $b^2\prod h_j\le T$, and for large $Y$ all $h_i-b$ with $i\ge3$ exceed $(m-3)(b+s)$. Together with $s>(m-2)a$ this gives (1).

There are $\gg_{m,\epsilon}(\log Y)^{m-3}$ choices of scales. The bounds $\delta_j H_j<h_j\le H_j$ allow each prime set at most $O_{m,\epsilon}(1)$ of them. Unique factorization recovers the ordered primes from $F$. Multiplication of the counts proves the stated bound.

## Verification

Run `python3 check_signed_power_pairs.py` with python-flint installed and `lean SignedDigitCarries.lean` with Lean 4.12.0. The Lean file checks the signed-digit step for arbitrary integer coefficient sequences satisfying the stated size and evaluation recurrences. The Python calculation checks the general criterion against binomial valuations, including a pair with bases $2$ and $65537$ at degree four.

For the general digit tests, 78,987 candidate cuts at small weights were checked by exact factorial-floor valuations. A separate dynamic program maximizes the number of positive parts avoiding both primes and checks (2) against that maximum. In the degree-four example all 386,850 proper candidate cuts are excluded.

An explicit five-prime case has $P=350708197319$ and offsets $(0,2,12,72,102)$. Taking $n=\prod_i(P+h_i)$, the prime $r=n-128$ completes the first three primes to a four-prime certificate at both $n$ and $n+1$. A row avoiding $r$ has two positive parts of total at most $129<P$. Such a row must be divisible in coefficient by all three of the base primes. The verifier certifies all six primalities. For this same product the three original primes certify two rows at **all** $dF,dF+1$ through $d=120663$ by (2). The uniform threshold covers $d\le102$ and the test using the scalar $L$ covers $d\le12563$. At $d=120664$ the bounds in (2) first fail for one of the three pairs.

For each of the 120,663 multipliers a fourth prime at most $dF$ gives a shared four-prime list for the two weights. The greatest distance to that prime is 1,962, at $d=45589$. The minimum sizes at the product and successor are $(4,4)$ for 118,612 multipliers, $(1,4)$ for 823, $(4,1)$ for 1,221 and $(1,1)$ for seven. The two minimum-size calculations use direct low-digit-sum tests and a separate enumeration by integer bisection. Large primes in this scan use the proved primality test in FLINT. [FLINT documentation](https://flintlib.org/doc/fmpz/fmpz_is_prime.html)

A different example keeps even the fourth prime fixed across a longer multiplier interval. [One four-prime list at $2^{76}$ weights](fixed_four_primes.md)
