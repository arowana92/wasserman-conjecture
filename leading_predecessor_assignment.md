# A large digit forces a row at an arbitrary cofactor

Let $p,q$ be distinct primes, $P=p^a$, $Q=q^m$ with $a,m\ge1$, $B$ a positive integer, $M=BQ$ and $N=1+PM$. The cofactor $B$ is unrestricted. Let $H$ be any power of $p$ with $H\le M$ and put $d=M-H$.

**If $d<2Q$, any counterexample must contain, up to order, the specific row**
$$
 \big(P(M-Q),\;Py,\;P(Q-y)+1\big),\qquad
          1\le y<Q,\quad Py\equiv1\bmod Q.                 \tag{1}
$$
This row must avoid both $p$ and $q$. Its assigned factor $U$ of the weight satisfies $3\le U\mid\gcd(B-1,y,N)$. In particular the full conjecture holds if $d<Q$, if this gcd is at most two, or if $(B-1)Q$ and $Q$ carry in base $p$. These criteria use no factorization of $B$ or $N$.

Use the [endpoint theorem](all_weights_endpoint_reduction.md) and assign every prime of $N(N-1)$ to a nonendpoint row avoiding it. Give the row taking $p$ the weight factor $U$ and the companion the factor $V$, so $UV=N$, $\gcd(U,V)=1$ and both factors are at least three. The parts of the first row are $P x,P y,P z+1$, where $x,y,z$ are positive and add without a base-$p$ carry to $M$. Indeed all the actual parts are multiples of $U$, excluding $z=0$. Because $M\ge H$, one slope receives a digit unit worth at least $H$. The other two slopes sum to at most $d$; at least one is unshifted and divisible by $U$, so $U<d$.

If the companion takes $q$, its three parts divided by $V$ sum to $U$ and two are positive multiples of $Q$, so $U>2Q$. This is impossible in the stated range. When the first row takes $q$, at least one of its two unshifted slopes is divisible by $Q$, by intersection of the zero positions. This slope is also divisible by $U$, prime to $Q$. It cannot be one of the two small slopes, since $UQ>d$. Its complement in $M$ is a positive multiple of $Q$ at most $d<2Q$, necessarily $Q$. The other unshifted slope now lies strictly between zero and $Q$, so its actual part must be one modulo $Q$. This gives (1). Dividing $U\mid (B-1)Q$ by the coprime factor gives $U\mid B-1$.

For example, if $B=1+p^i q^j$ with $i,j\ge0$ then the gcd condition always holds. With $B=55$, $p=2$, $q=3$ this covers every $a\ge1$ whenever $53\cdot3^m<2^n\le55\cdot3^m$ for some $n$. There are infinitely many such $m$, of natural density $\log_2(55/53)$, by uniform distribution of the fractional parts of $m\log_2 3$.

## A larger window in which both primes must go to the same row

At every $B$, **a split assignment requires $d>5Q$**. In fact write two actual parts on the row taking $p$ as $UPi,UPj$ and two on the row taking $q$ as $VQk,VQl$. Here $i,j,k,l$ are positive. The positive remaining parts give $V>P(i+j)$ and $U>Q(k+l)$, whence $(i+j)(k+l)<B$ by integrality and $UV=BPQ+1$. More explicitly $UV\ge(P(i+j)+1)(Q(k+l)+1)>PQ(i+j)(k+l)+1$.

Also $B\mid ijkl$. Each full power in $B$ supported away from $p,q$ divides a displayed multiplier on the row to which its prime is assigned, by intersection of the two zero positions. The $p$-part of $B$ divides both $i,j$ because the two unshifted parts are divisible by the full $p$-power in $N-1$. The $q$-part similarly divides $k,l$. [Predecessor assignments at arbitrary cofactors](predecessor_digit_splits.md)

At least one of $Ui,Uj$ is among the two small slopes, so $d>Q\min(i,j)(k+l)$. If $\min(i,j)\ge3$ the integer multiplier on the right is at least six. If the minimum is two and $k+l=2$, then $(i+j)(k+l)>ijkl$, impossible. If the minimum is one and $k+l\le4$, then $kl\le k+l$ and again the same contradiction follows. This proves the bound.

Consequently for **each fixed $B\ge6$ and distinct primes $p,q$**, a set of exponents $m$ of lower natural density at least
$$
                 \min\big(1,\log_p(B/(B-5))\big)
$$
certifies the conjecture at **every** $1+B p^a q^m$, $a\ge1$. Indeed the range $(B-5)Q\le H\le BQ$ has the indicated density by irrational rotation. With both primes on one row some unshifted slope equals $tQ$, $1\le t<B$. Thus $t q^m$ and $(B-t)q^m$ must add without a base-$p$ carry. For each fixed $t$ this permits only $O_{B,p}(X^{\beta_p})$ indices $m\le X$, where $\beta_p<1$, by the [exponent digit bound](three_factor_exponent_rectangles.md#how-many-exponent-terms-can-pass-a-cut). If $p\mid B$, first divide the slopes $t,B-t$ by their common $p$-power. If both reduced slopes are units and their sum is zero modulo $p$, they carry at the first occupied place for every $m$. In the other cases the same digit-block proof applies. Summing over the $B-1$ cuts removes a set of density zero.

At $B=35$ the split multipliers and the low and high binary digits give a much stronger lower density for $p=2$. [Cofactor 35](cofactor_35.md)

Run **python3 check_leading_predecessor_assignment.py**. Direct enumeration at $M\le240$, five prime bases and five shifts per base finds 47,774 assignments passing the first size and digit conditions. Twelve have the required residues modulo $Q$. All twelve equal (1), and four actually avoid both primes, as checked by exact binomial products. For example, at $N=573$ the forced row $(12,528,33)$ avoids 2 and 11. With the endpoint and the companion $(191,191,191)$ its coefficient gcd is 39.
