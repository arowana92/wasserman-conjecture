# A square-root window for two predecessor powers

Let $P=p^a<Q=q^b$ be powers of distinct primes with $a,b\ge1$, and write $h=Q-P$. For every integer $D\ge1$,
$$
 8P>D^3h^2+4D^2
 \quad\Longrightarrow\quad
 \text{Wasserman's three-row conjecture holds at }N=DPQ+1 .
$$
The multiplier may share either base. At fixed $D$ this gives gaps of order $\sqrt P$. The earlier [general cluster theorem](cluster_factor_bounds.md) at two powers gave gaps of order $P^{1/3}$.

Use the [endpoint reduction](all_weights_endpoint_reduction.md) and suppose a failure. The other two rows have assigned weight scales $U,V\ge3$ with $UV=N$, each dividing all parts of its row. Choose $A$ to take $p$, with scale $U$. At $D=1$ the [two-factor result](endpoint_reduction.md) already applies, so take $D\ge2$. Our bound gives $Dh<P$ and $D<P$. Grouping the carry-free base-$p$ digits of $A$ into blocks of size $P$ gives three expressions
$$
 A_i=a_iP^2+b_iP+\delta_i,\qquad \sum a_i=D,\quad \sum b_i=Dh,\quad \sum\delta_i=1
$$
with nonnegative coefficients. Here $N=DP^2+DhP+1$.

## When the primes go to different rows

If the companion takes $q$, two of its normalized parts are positive multiples of $Q$. Their sum is strictly less than $U$, so $U>2Q$. Of the two parts of $A$ with $\delta_i=0$, select the one with smaller middle coefficient and denote its coefficients by $a,b$. Then $b\le Dh/2$. We have $U\mid aP+b$ since $P$ is a unit modulo $U$. Reducing $a^2N$ modulo $U$ gives
$$
 U\mid a^2+D b(b-ah),\qquad |a^2+D b(b-ah)|\le D^2+D^3 h^2/4<2Q .
$$
For the upper bound use $b^2\le D^2h^2/4$. For the lower bound complete the square to get $b(b-ah)\ge-a^2h^2/4$. Consequently the displayed integer is zero. Neither $a$ nor $b$ can vanish, since $A_i>0$. Write $a=gk,b=gl$ with $(k,l)=1$. Then $Dl(kh-l)=k^2$ with both factors $l,kh-l$ positive and coprime to $k$. Both equal one. Hence $kh=2$ and $D=k^2$. The only possibility with $D\ge2$ is $(D,h)=(4,1)$. In this case one of the consecutive powers is even, so $4PQ$ has at most two distinct prime factors and the two-factor theorem applies.

## When one row takes both primes

At least one part is divisible by both powers, by intersecting the two sets of zero residues. Its value is $xPQ$ with $0<x<D$, since the other two positive parts have sum at least two. Group the remaining parts, then remove their units digit in base $q$. Dividing by $Q$ leaves the carry-free cut $xP,(D-x)P$. This is the [binomial transfer](clustered_predecessors.md#a-binomial-test-and-passage-to-the-next-weight). In fact for every $q$-power $L$ and $s=[-P]_L$ (least nonnegative residue), it requires
$$
 Ds>L. \tag{1}
$$
Indeed $s>0$. If $Ds\le L$, the two residues are $L-xs,L-(D-x)s$, which sum to at least $L$. Taking $L=Q$ gives $s=h$ and contradicts $Dh\le Q$. This finishes the proof of the square-root window.

**Linear window at certain cofactors.** For $D$ a prime power, twice an odd prime power, or $1,12,15$, the full conjecture holds whenever $D(Q-P)\le Q$. More generally any violation of (1), in either direction, suffices. After absorbing into $P,Q$ all factors of $D$ sharing their bases, the remaining cofactor is a divisor of $D$; every such divisor also belongs to this list. The [four-bin theorem](two_prime_powers.md) therefore forbids splitting the two bases. The same conclusion holds at $D=21,33$ by the [parity refinement](four_factor_squares.md#why-the-two-bases-must-be-assigned-to-the-same-row). The common-row proof uses the original powers and original $D$, whether or not factors were absorbed. For example at $D=10$ all pairs with $9Q\le10P$ are covered.

## Counting possible exceptions throughout wider windows

For each fixed $D$ and $\eta>0$, there are at most $O_{D,\eta}(h^{1+\eta})$ possible exceptional pairs at any fixed gap $h\ge1$, even if $P$ is unbounded. Thus the total through gap $H$ is $O_{D,\eta}(H^{2+\eta})$.

Indeed $P\le Dh$ accounts for at most $Dh$ pairs. For $P>Dh$, (1) forbids grouping. Take the split row at $P$ and a part with units block zero as above. Here $0\le a\le D$, $0\le b\le Dh$, and $U>2Q$. Necessarily $a\ge1$, since otherwise the positive $b$ would be divisible by $U$ despite $b<U$. Fix $h,a,b$. The nonzero integer $a^2+D b(b-ah)$ has absolute value $O_D(h^2)$, hence only $O_{D,\eta}(h^\eta)$ positive divisors $U$ by the elementary divisor bound. The zero cases have already been excluded. For each divisor, $aP+b=iU$ with $1\le i\le D$ since $aP+b\le(D+1)P$ and $U>2P$. Hence only $O_D(1)$ choices of $P$ remain. Summing over $a,b$ proves the bound.

In particular fix $0<\theta<1$ and take $H=\lfloor Y^\theta\rfloor$. Among pairs of ordinary primes in $[Y,2Y]$ at distance at most $H$, the proportion of possible exceptions tends to zero. To see the denominator, partition the interval into $O(Y/H)$ bins of length at most $H$. The prime number theorem and Cauchy–Schwarz give $\gg YH/(\log Y)^2$ distinct pairs within bins, and one can choose $\eta>0$ with $\theta(1+\eta)<1$ in the numerator.

## Three nearby powers

The same normalized-size observation improves the successor bound for three powers in the [general cluster theorem](cluster_factor_bounds.md). Suppose $P,P+s,P+h$ are powers of distinct primes, $0<s<h$, and put $H=h+1$. The conjecture holds at $D P(P+s)(P+h)+1$ whenever
$$
 P>(12D)^{20}H^3 .
$$
Indeed in that proof take $n=3$ and let $k$ be the number of the other powers assigned to the companion. Its coefficient and translation bounds only require $P>4(D(2H)^3+1)+4$, implied here. When $k=0$ the negative coefficients at the largest translated node give the contradiction there. When $k=1,2$, at least one positive normalized part of the companion is divisible by the product of *all* its $k$ powers: the zero pairs of two powers intersect in a triple. Hence $U>P^k$. The resultant bound in the same proof is $U\le K H^{3k}$ with $K=(2\cdot3^3)^3D^6$, and the cubic weight polynomial cannot be a square. Our threshold exceeds $K H^3$, a contradiction. This improves the general gap exponent from $2/9$ to $1/3$ at three powers. Separate bounds for particular gap orderings and multipliers are in the [three-factor note](clustered_predecessors.md).

Run **python3 check_predecessor_pair_gaps.py** for the integer identities, their zero cases, and direct carry-free row comparisons. [Verifier](check_predecessor_pair_gaps.py)
