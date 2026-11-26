# Three prime powers in a short interval

Let $P<Q<R$ be powers greater than one of three distinct primes $p,q,r$. Write $M=PQR$, $s=Q-P$, $t=R-Q$ and $h=R-P$.

**Two consecutive weights at every multiplier.** Suppose $t>s$. For every positive integer $d$ with $2dh^2\le P$, Wasserman's trinomial conjecture holds at **both $dM$ and $dM+1$**. At each of those weights any two positive trinomial coefficients already share a prime from $\{p,q,r\}$. At $dM$ this last assertion even holds for any two nontrivial binomial coefficients. The multiplier may share primes with $M$.

**Successors for either ordering of the gaps.** For $2h^2\le P$, the full conjecture holds at $M+1$, without the assumption $t>s$. If all three powers are odd and $4h^2\le P$, it also holds at $2M+1$ for either ordering.

The first result gives $\gg_\epsilon X^{2/3}/(D\log^3 X)$ distinct products $M\le X$ of three primes at which *every* $1\le d\le D$ works simultaneously, uniformly for integers $1\le D\le X^{1/3-\epsilon}$. The proofs below use digit comparisons at two powers of each base and the [endpoint theorem](all_weights_endpoint_reduction.md).

## A binomial test and passage to the next weight

Say a row avoids a prime when its coefficient is coprime to that prime. Let $W<Z$ be powers of distinct primes dividing $n$, and take powers $L_W,L_Z$ of their respective bases. Suppose there is a positive multiple $K$ of $Z$ with
$$L_W\le WZK<L_Z.$$
Set $D=\lfloor n/L_W\rfloor$, $E=\lfloor n/L_Z\rfloor$, $H=n-DL_W$, $m=n/(WZ)$, and $\delta=[-W]_Z$ (least nonnegative residue). No positive binomial row at $n$ can avoid both base primes if either of the following holds:
$$D=E,\qquad \delta(m-KD)\le Z; \tag{1}$$
$$E\ge D-1,\qquad \delta\lfloor H/(WZ)\rfloor+[n/Z]_Z<Z. \tag{2}$$

Indeed, by carry freedom each part is $A_i=WZ f_i$ with $f_i>0$. Put $u_i=\lfloor A_i/L_W\rfloor$ and $e_i=f_i-Ku_i$. At each prime-power layer, the floors of the parts add to the floor of the weight. Consequently $\sum u_i=D$, the remainders at $L_W$ sum to $H$, and every positive $e_i$ is at most $H/(WZ)$.

In (1) the floors at $L_Z$ are termwise at most those at $L_W$ and have equal total, so both are $u_i$ at position $i$. Thus $e_i>0$ at each position, also when $u_i=0$ by positivity. Dividing the parts by $Z$ preserves carry freedom in its base. Their low residues at $Z$ become $[-\delta e_i]_Z\ge Z-\delta e_i$. Their sum is at least $2Z-\delta(m-KD)\ge Z$, a contradiction.

In (2), any positive $e_i$ would give a residue strictly larger than the entire budget $[n/Z]_Z$. Thus both $e_i\le0$. Then $u_i\ge1$ and $A_i\le WZKu_i<L_Zu_i$, so their two floors at $L_Z$ sum to at most $D-2$, also a contradiction.

**Transfer.** When $WZ\mid n$ and the two primes are incompatible for positive binomial rows at $n$, they are incompatible for positive trinomial rows at **both** $n$ and $n+1$. At $n$, group two parts of any simultaneously avoiding row. At $n+1$ the low residues at each of $W,Z$ are $0,0,1$. Some part $a$ is divisible by both. Group the others and subtract one from their sum. The resulting pair $(a,n-a)$ is positive. In either base the grouped second part was 1 modulo the corresponding power, so removing its units digit preserves carry freedom. This gives the forbidden binomial row at $n$.

Three pairwise incompatible primes divide the endpoint coefficient $N(N-1)$ at each of $N=n,n+1$. Two rows can avoid at most two of these primes. The endpoint theorem then excludes every three-row counterexample at either weight.

## Checking the three pairs

Take $n=dM$, $2dh^2\le P$, and $t>s$. We have $PR-Q^2=(t-s)Q-st>0$, so $P^3<Q^3<M<R^3$. Both $n/P^3$ and $n/Q^3$ have integer part $d$, and $n/R^3$ has integer part $d-1$. In fact
$$ d(M-P^3)/P^3\le d(2h/P+h^2/P^2)<1,\qquad
 d(R^3-M)/R^3<3dh/R<1,$$
using $h\ge2$. For each ordered pair of powers $W<Z$ take $L_W=W^3$, $L_Z=Z^3$, $K=Z$, and $\delta=Z-W$.

For $W=P,Z=Q$, condition (1) follows from $\delta(m-Kd)=dst<Q$. For either pair with $Z=R$, write $T$ for the remaining power. We have $[n/Z]_Z=d(Z-W)(Z-T)$ since this product is less than $Z$. With $H=n-dW^3$, the left side in (2) is at most
$$d(Z-W)\big(T-W^2/Z+Z-T\big)
 =d(Z-W)^2(1+W/Z)<2dh^2\le P.$$
This proves the first theorem.

The gap ordering matters for this particular set of three primes. For example, when $P,Q,R$ are odd primes, $2h^2\le P$ and $3\le t\le s$, take $f=\lfloor R/t\rfloor$. The row $QR(f,f,P-2f)$ at $M$ avoids **both** $Q$ and $R$. To verify this, $P-2f>0$ since $t\ge3$ and $R<3P/2$. Put $b=R-tf$, so $1\le b<t$. After division by $Q$, the two identical entries $Rf$ have low digit $t-b$ at $Q$, within the total budget $[PR]_Q=Q-st$ even when taken twice. After division by $R$, the corresponding two low digits are $b,b$ with total at most $[PQ]_R=ht$. Both normalized totals are less than their respective base squares since $t\le s$. Passing the low digits suffices for the whole addition. Each individual prime $P$ is avoidable here because the odd composite weight has odd base-$P$ digit sum and is not a power of $P$.

## Rounding the CRT multipliers at a successor

There is a general integral restriction behind the result at $M+1$. In a possible endpoint counterexample write $UV=N$, $ST=N-1$ for the products assigned to rows $A,C$ from the weight and the predecessor, respectively. Partition $S=w_1w_2w_3$ by the exceptional position of each assigned prime power. Thus $(S/w_i)\mid A_i$, $U\mid A_i$, and $\gcd(U,N-1)=1$. Pick any positive integers $m_i$ divisible by $w_i$. If $U\nmid m_j$ and $A_j\ge N-H$, then
$$S\le m_j+H\max_i|m_i-m_j|. \tag{3}$$
Indeed $B=\sum m_i A_i/S$ is an integer divisible by $U$. Hence $B-Tm_j$ is a nonzero integer. Multiplying by $S$ gives $m_j+\sum_{i\ne j}(m_i-m_j)A_i$, proving (3).

For the particular successors $N=dM+1$ under consideration, use $d=1$ or $d=2$. When $d=1$, the [three-factor assignment](endpoint_factors.md) gives a row $A$ avoiding every base of $M$, with those powers exceptional at distinct positions. For $d=2$ with $P,Q,R$ odd, the [four-factor assignment](two_prime_powers.md) forces these three odd bases onto one row. Their exceptional positions are distinct since a position exceptional for none would be divisible by $M$ and by a coprime scale $U\ge3$, exceeding $2M+1$. In both cases write $w_i=P,Q,R$ in the appropriate order and $A_i=(M/w_i)k_i$, where $U$ divides all three integers $k_i$.

Let $b_i=\lfloor A_i/P^3\rfloor$ and $z_i=A_i-b_iP^3$. Under $2dh^2\le P$ the weight is between $dP^3$ and $(d+1)P^3$; this follows from the preceding estimate with the extra $1/P^3$. Carry freedom at the base of $P$ gives $\sum b_i=d$, $\sum z_i=H=d(M-P^3)+1$. Then
$$M\left(\sum k_i-\sum b_i w_i\right)
 =P+\sum(w_i-P)z_i-(M-P^3)\sum(w_i-P)b_i.$$
Its absolute value is at most $P+hH<M$. For the strict bound use $(M-P^3)/M<(s+h)/P$, $dh(s+h)\le P-dht$, and $P+h=R<(dht/P)M$. Thus the sum in parentheses vanishes.

For $d=1$ the right sum $\sum b_i w_i$ is one of the three powers and is coprime to $U$, while the left sum is a multiple of $U\ge3$. For $d=2$ all three $w_i$ are odd. Now $k_i$ and $A_i$ have the same parity, so $\sum k_i$ is odd and $\sum b_i w_i$ is even. These contradictions prove the two successor statements.

## Three powers on different scales

The same rounding bound also works with widely separated factors. Let $B,X,Y$ be positive powers of one prime $p$, let $a,b>0$ be integers, and suppose $X+a$ and $Y+b$ are prime powers. Set $C=\max(B,X,Y)$ and $\eta=\max(a/X,b/Y)$. The full conjecture holds at
$$N=1+B(X+a)(Y+b)$$
under either of these sufficient conditions:
$$8C\eta^2\le1,\qquad\text{or}\qquad
 B\le\min(X,Y)\ \text{ and }\ 8ab\eta\le1.$$

When the predecessor has at most two distinct prime bases the [two-factor result](endpoint_reduction.md) applies. Otherwise put $M=N-1$ and $L=BXY$, and take the mandatory CRT row. The relative gap satisfies $(M-L)/M<2\eta$. In either case $\eta\le1/4$, so $L>N/2$ since $L\ge8$. The row avoiding $p$ has one part at least $L$. Under the first condition use in (3) the three multiples $C,\ C(X+a)/X,\ C(Y+b)/Y$. Each multiplier of its corresponding predecessor factor is a power of $p$ or one, hence all are coprime to $U$. Here $S=M$, the largest multiple is at most $C(1+\eta)$ and the spread at most $C\eta$. The upper bound in (3), divided by $M$, is strictly less than $C(1+2\eta)/M+2C\eta^2<1$ since $M>4C$.

Under the second condition, $B$ divides $X,Y$. The CRT multiplier $k$ at $B$ satisfies $kab\equiv1\bmod B$, $1\le k<B$, $U\mid k$. If $ab=1$, this already makes $k=1$. Otherwise $ab\min(k,B-k)\ge B-1$ since neither signed product with a positive multiplier can equal $1$ in absolute value. The exceptional part is $Mk/B$. If it received the leading $L$ then $M-Mk/B\le N-L$; if another part received it then $Mk/B\le N-L$. Thus $(B-1)/B\le ab(N-L)/M<2ab\eta+ab/M$. But $2ab\eta\le1/4$ and $ab/M<\eta^2/B$, contradicting $B\ge2$.

For example the first condition works at $(B,X,Y)=(2^{32},2^{24},2^{32})$ with $(a,b)=(43,15)$. The second works at $(2^{20},2^{20},2^{90})$ with $(a,b)=(7,133)$. In both cases the shifted factors are prime.

## Counting whole intervals of multipliers

Take primes in $[Y,2Y]$ for $Y=X^{1/3}/2$. Divide into bins of length at most $H=\lfloor\sqrt{Y/(2D)}\rfloor$ using $O(Y/H)$ bins. Under $D\le X^{1/3-\epsilon}$ the length $H$ grows as a positive power of $Y$. By the prime number theorem, bins with at least $cH/\log Y$ primes account for $\gg Y/\log Y$ primes, for a sufficiently small fixed $c>0$.

In one such bin write $k$ for its prime count. An interval of length $\lambda H$ contains at most $2\lambda H/\log(\lambda H)+1<k/4$ primes for a sufficiently small fixed $\lambda>0$ depending on $\epsilon$. This uses the short-interval Brun–Titchmarsh inequality. Select $P,Q$ from the bottom quarter of the ordered primes in this bin, in the same subinterval of length at most $\lambda H$. Partitioning the bin into $O_\epsilon(1)$ such subintervals gives $\gg_\epsilon k^2$ choices. Any $R$ in its top quarter satisfies $R-Q>\lambda H\ge Q-P$, since there are at least $k/2-O(1)$ primes between the two quarters. There are therefore $\gg_\epsilon k^3$ suitable triples per bin. [Montgomery and Vaughan, equation (1.12)](https://personal.science.psu.edu/rcv4/personal/Publications/large_sieve.pdf)

Summing and using $k^2\gg H^2/\log^2Y$ in the counted bins gives $\gg_\epsilon YH^2/\log^3Y$ triples. Unique factorization gives distinct weights $M\le X$, and $2Dh^2\le Y$ in each bin. This proves the quantitative assertion.

The [shared four-prime certificates](shared_cluster_certificates.md) strengthen this count for smaller growing multiplier ranges and give the exact verification commands.
