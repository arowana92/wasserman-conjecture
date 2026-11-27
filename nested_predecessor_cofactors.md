# Four cofactors with unrestricted nested predecessor powers

**Theorem.** Let $P<Q$ be powers of distinct primes with positive exponents and $P\mid Q-1$. Wasserman's full three-row conjecture holds at
$$
                 N=1+DPQ,\qquad D=5,8,9,10 .
$$
There are no exponent or size restrictions. At $D=10$ this includes all the square and nonsquare outer powers left by the [earlier nested-power results](ten_nested_powers.md). For each fixed $P$, Dirichlet's theorem supplies infinitely many ordinary primes $Q$ in the indicated progression.

At $D=7$ a possible counterexample requires $P\equiv55\bmod56$, $Q\equiv19\bmod56$ and the two nonendpoint rows
$$
 (4PQ,(2P-1)Q+1,(P+1)Q),\qquad (N/4,N/4,N/2).
$$
In particular a square on either side or a cube on the $Q$ side certifies this cofactor as well.

## The assigned row

Write $W=PQ$, and take $b=D$ except that $b=5$ at $D=10$. The [static numerator theorem](static_square_cofactors.md) and [distinct exceptional positions](predecessor_inverse_digit_gaps.md) imply that in a counterexample $(b,W)=1$. One nonendpoint row $A$ must avoid the bases of $b,P,Q$. A unitary divisor $U\ge3$ of $N$ divides all its parts, and $A$ avoids every prime of $U$. Its part exceptional at $b$ is $xW$, where $U\mid x$, $x<D$, $(x,b)=1$. These results apply also when a varying power at cofactor ten is even. They use the [all-weight endpoint reduction](all_weights_endpoint_reduction.md).

More generally, let $s$ be the least positive inverse of $Q$ modulo $P$, and put $E=sQ$. The part exceptional at $P$ is a multiple of $Q$ congruent to one modulo $P$. The row therefore has the form
$$
 A=(xW,\ tW-E+1,\ kW+E),\qquad x+t+k=D,\quad k\ge0,\quad t\ge1 .
$$
The last part is divisible by $bU$, hence $bU\mid kP+s$. When $P\mid Q-1$, we have $s=1$, $k\ge1$ and $(k,bU)=1$. If $U<7$ then $A$ must also avoid 2: the three normalized parts of the companion sum to $U$ and cannot have distinct 2-adic valuations unless their sum is at least $1+2+4$. The endpoint coefficient is always even.

## A finite leading-digit test at unbounded weights

Suppose $C\ge1$ and $1\le E\le W/C$. For a required prime $r$ let $H$ be the largest power of $r$ at most $W$ (allowing $H=1$), and put $w=W/H\in[1,r)$. If for any $j\ge0$ one has
$$
 \left\lfloor\frac{xw}{r^j}\right\rfloor+
 \left\lfloor\frac{tw}{r^j}\right\rfloor+
 \left\lfloor\frac{(k+1/C)w}{r^j}\right\rfloor
 <
 \left\lfloor\frac{Dw}{r^j}\right\rfloor , \tag{1}
$$
then the row carries at $r$. The left side bounds the sum of the parts' floors at $Hr^j$ from above and the right side bounds the floor of the total from below. A strict deficit contributes positively in Legendre's formula. At modulus one a strict deficit cannot occur.

The table lists all the nested possibilities and a prime and a value of $C$ giving (1) throughout $[1,r)$ for some $j\in\{0,1,2,3\}$, except for the three entries treated separately. Enumerate $3\le U\mid x$, $(U,D)=(x,b)=(k,bU)=1$ and $x+k<D$ to obtain the list.

| $D$ | $(U,x,t,k)$ | $(r,C)$ in the same order |
|---:|---|---|
| 5 | $(3,3,1,1)$ | $(2,4)$ |
| 7 | $(3,3,3,1),(3,3,2,2),(5,5,1,1)$ | $(2,1),(2,1),(2,6)$ |
| 7 | $(4,4,2,1)$ | Residual row |
| 8 | $(3,3,4,1),(5,5,2,1)$ | $(2,1),(2,1)$ |
| 9 | $(4,4,4,1),(5,5,3,1),(5,5,2,2)$ | $(2,1),(3,2),(2,4)$ |
| 9 | $(7,7,1,1)$ | Repunit argument |
| 10 | $(3,3,6,1),(3,3,5,2),(3,3,3,4)$ | $(3,2),(3,1),(3,1)$ |
| 10 | $(3,6,3,1),(3,6,2,2),(7,7,2,1)$ | $(3,2),(2,2),(5,9)$ |
| 10 | $(7,7,1,2)$ | Repunit argument and $P=17$ |

In every leading-digit entry $bU\mid kP+1$ rules out $2\le P<C$. Hence $E=Q\le W/C$. To see the small size of this computation, consider the remaining square-$Q$ row at ten, $(x,t,k)=(6,3,1)$. At $r=3,C=2$, a point escaping (1) at $j=3$ must lie in $[1,27/10)$. At $j=2$ this shrinks to $[3/2,9/5)$ and none of that interval escapes at $j=1$.

For all 15 entries there are two exact verifications. The first partitions $[1,r)$ at the discontinuities of the four floors for all four layers, then checks their values at the left endpoints. The second takes $M=\operatorname{lcm}(x,t,Ck+1,D)$ and $z=\lfloor Mw\rfloor$. Each discontinuity lies on the grid of spacing $1/M$ because $r^j$ is integral. It suffices to check each integer $M\le z<rM$ at $w=z/M$. In total the first calculation uses 294 rational intervals and the second uses 964 grid points. Lean checks the grid reduction for integer ratios and all 15 finite covers. [Python checks](check_nested_leading_layers.py) · [Lean certificate](NestedLeadingLayers.lean)

## The repunit case

For any $D$, a nested row of the displayed form with $t=1$ and $P>Dk$ cannot have a weight scale $U>1$. Write $Q=\sum_{j\ge0}v_jP^j$ in base $P$, with $v_0=1$, and $m=(Q-1)/P$. Above the $P$ shift the three parts are $xQ,Q-m,kQ+m$, adding without a carry in the base prime of $P$. Thus they also add without carries across or within the sequence of base-$P$ blocks. At position $j$ the raw coefficients are
$$
              x v_j,\qquad v_j-v_{j+1},\qquad k v_j+v_{j+1}.
$$
Assume the digits through $v_j$ are all one. Earlier raw coefficients are ordinary digits, so no multiplication borrow or carry enters this position. If $1<v_{j+1}<P-k$, the negative middle coefficient adds an extra $P$ upon reduction, pushing the digit sum above its budget. If $v_{j+1}=P-l$ with $1\le l\le k$, group the last two upper parts. At position $j+1$ the two numbers $xQ,(D-x)Q$ have residues $P-xl,P-(D-x)l$. Lower digits were all one, and $Dl<P$. The residues therefore force a carry. Consequently the next digit is zero or one. After a first zero, any next nonzero $v$ gives raw digits $0,-v,v$ and forces a carry. The digits hence form a repunit $Q=1+P+\cdots+P^n$. The second actual part is $(P-1)Q+1=P^{n+1}$, impossible because $U>1$ divides it and $(U,P)=1$.

At nine, $63\mid P+1$ implies $P>9$. At ten, $35\mid2P+1$ leaves only $P=17$ at $P\le20$. Here $Q\equiv4\bmod5$, since $7PQ\equiv1\bmod5$. Write $Q=4+5v$. In the parts $119Q,16Q+1,35Q$ the second base-five digits are $4v,3+v,3$ reduced modulo five. Their sum is six or eleven. This proves the four unrestricted assertions.

At seven the remaining scale is four, so the positive normalized companion is $(1,1,2)$. We have $28\mid P+1$ and $4PQ\equiv1\bmod7$, hence $Q\equiv5\bmod7$. Unitarity gives $4\parallel N$, so $PQ\equiv5\bmod8$. Both powers are odd: $P\equiv3\bmod4$ and $PQ\equiv5\bmod8$. If $P\equiv3\bmod8$, the parts $4PQ$ and $(P+1)Q$ have the same lowest binary one, impossible since $A$ avoids 2. Thus $P\equiv7\bmod8$ and $Q\equiv3\bmod8$, proving the stated residues.

## Inverse residues beyond nesting

At cofactor ten let $P,Q$ be powers of distinct primes in either order and $s=[Q^{-1}]_P$. If $(s,15)=1$, the full conjecture also holds under either of these conditions:
$$
 2s\le P,\quad PQ\not\equiv23\bmod35;
 \qquad\text{or}\qquad
 9s\le P,\quad 35\nmid 2P+s .
$$
Indeed $5\nmid s$ forces $k>0$ and $5\nmid k$. At $U=3$ also $3\nmid k$, giving exactly the five leading-digit entries for this scale. All need only $C\le2$. At $U=7$ one has $x=7$, $k=1,2$, and $PQ\equiv23\bmod35$ by $7\mid N$ and $7PQ\equiv1\bmod5$. If $9s\le P$, (1) excludes $k=1$, leaving $35\mid2P+s$ for $k=2$. These arguments use the original powers even when one is even.

A further test applies when $D=b$ is any prime power or $D=2b$ with $b$ an odd prime power. Write $P=p^a$. For any layer $L=p^j$, $j\ge1$, put $z=[DQ]_L$. In a possible failure with $Dz<L$, necessarily $b\mid z$. In fact $xQ,(D-x)Q$ add without a carry after the shift by $P$. Hence $v=[xQ]_L\le z$. Since $Dv\equiv xz\bmod L$ and both sides are in $[0,L)$, $Dv=xz$. Use $b\mid D$ and $(b,x)=1$.

For example at ten the conjecture holds also whenever $P\mid sQ-1$ for $s=2,5,10$. For $P>10$ and $s\mid10$, the upper cut at $L=P$ gives $v=[xQ]_P\le 10/s$. Thus $sv\equiv x\bmod P$ and $0\le sv\le10<P$, giving $sv=x$. At $s=5,10$ this contradicts $5\nmid x$. At $s=2$ it forces $x=6$, and the cofactor-ten leading tests for $U=3$ apply since $2s\le P$. For $P\le10$ write $N=1+dq^e$ using the full varying $q$-power in $10PQ$. The remaining cofactor $d\le100$ is coprime to $q$, so the [verified small-cofactor theorem](prime_power_cofactors.md) and the all-weight endpoint reduction apply.

Reproduce the finite arithmetic with **python3 check_nested_leading_layers.py** and **lean NestedLeadingLayers.lean** using Lean 4.36.
