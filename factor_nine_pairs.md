# Two of the three weights $5t+1$, $7t+1$, $9t+1$

**Let $t=p^i q^j$ for primes $p,q$ and integers $i,j\ge1$. At least two of $5t+1$, $7t+1$ and $9t+1$ satisfy Wasserman's full trinomial conjecture.** The statement allows $p=q$.

For even $t$ the proof gives a test to identify the guaranteed weights. Put $u=t/2^{v_2(t)}$. Among the $315\cdot4=1260$ pairs $(t\bmod315,u\bmod8)$, only 21 can admit a factor-five exception, 63 a factor-seven exception, and 19 a factor-nine exception. These three sets are disjoint. The other 1,157 residue pairs certify all three weights.

**For the factor-nine family alone**, in any interval of $H\ge2$ consecutive weights the number of possible exceptions is $O(H^{0.12})$. The bound improves to $O_\epsilon(H^\epsilon)$ for every $\epsilon>0$ outside the single case
$$
 t\equiv115\bmod144,\qquad
 V=(9t+1)/4\ \text{has no prime factor at most }100{,}000 .
$$
In particular the stronger bound holds for all the odd weights in this family.

## The allowed rows at factor nine

If $9t$ has fewer than three distinct prime factors, the [two-factor theorem](all_weights_endpoint_reduction.md#consequences) applies. Otherwise 9 is the full power of 3 and the [three-factor assignment](endpoint_factors.md#when-n-1-has-three-distinct-prime-factors) gives a nonendpoint row $A$ taking all three predecessor primes. It has a part $at$, with $3\le U\mid a<9$, $at\equiv1\bmod9$, $UV=9t+1$, and $\gcd(U,V)=1$. Here $U$ is the product of the weight powers assigned to $A$. The other nonendpoint row $C$ has parts $V y_i$ with the three positive $y_i$ summing to $U$. Thus $(a,U)$ is one of $(4,4),(5,5),(7,7),(8,4),(8,8)$.

When $t$ is even the weight is odd and only $(5,5),(7,7)$ remain. For an even input at factor $w=5,7,9$, each possible shape has a cut $at,(w-a)t+1$ that must be binary carry-free. The units digit 1 lies below the shift $2^{v_2(t)}$, so $au$ and $(w-a)u$ have disjoint bits. Their lowest three bits suffice for the table. Retain a residue pair if any listed $(a,U)$ satisfies $at\equiv1\bmod w$, $wt+1\equiv0\bmod U$, and the binary test. The factor-five and factor-seven shapes follow from the [two-family argument](five_seven_pairs.md).

| $w$ | $(a,U)$ | Allowed residue pairs |
|---:|---|---:|
| 5 | $(3,3)$ | 21 |
| 7 | $(3,3),(6,3),(5,5)$ | 63 |
| 9 | $(5,5),(7,7)$ | 19 |

All moduli for $t$ divide 315. Python counts each set by explicit bit intersection, and tests their pairwise intersections. Independently, Lean checks that at each residue pair at most one family can pass. Its binary test requires no overflow at moduli 4 and 8. The units cannot overlap since $w$ is odd and one summand is even. This proves the assertion at every even input. [Kernel-checked residue table](CarryOrder.lean)

For odd $t$, $5t+1$ satisfies the conjecture by the [even-weight theorem for factor five](predecessor_five_seven.md#full-factor-five). If $7t+1$ is exceptional, the [even factor-seven restriction](predecessor_five_seven.md#a-restriction-on-every-base-seven-digit) gives $t\equiv11\bmod72$. Then $t\equiv2\bmod9$ forces $a=5$, $U=5$ at any factor-nine exception. Since $9t+1$ is even, $C$ must avoid 2. The $y_i$ consequently have distinct 2-adic valuations, because two equal valuations give an overlapping bit after scaling by $V$. Their sum is at least $1+2+4=7$, contradicting $U=5$. This completes the proof for odd inputs.

## Counting factor-nine exceptions in moving intervals

Use the [uniform intersection bound for digit sets](five_seven_pairs.md#from-the-digit-sets-to-moving-intervals). For $t$ even, $A$ must avoid $2,3,U$. At block length 256 the dimensions of the digit sets for $at,(9-a)t+1$ are bounded in the three respective bases by $(0.59,0.56,0.71)$ when $a=5$ and $(0.58,0.50,0.77)$ when $a=7$. Both sums are below 2.

For odd $t$ consider the following possibilities.

* If $a=8$, then $2\mid U$ and the same cut must avoid 2 and 3. The two block dimension bounds are $0.71$ and $0.07$.
* If $a=5$, the argument with $y_1+y_2+y_3=5$ rules out the shape. If $a=7$, $C$ must avoid 2, and its normalized parts of sum 7 must be $1,2,4$ by the distinct valuations. The binary graph for $V,2V,4V$ has a dimension bound of $0.56$. The ternary and base-seven cut bounds for $t$ are $0.50$ and $0.77$. In the intersection argument express $t=(7V-1)/9$ as an affine function of $V$, which only multiplies the bounded binary slope by $7/9$ and changes the translation.
* If $a=4$, then $U=4$ is unitary, so $4\parallel9t+1$ and $t\equiv3\bmod8$. The cut $4t,5t+1$ avoids 2. The two bits of $4t$ at positions 2 and 3 force $5t+1$ to be 0 modulo 16, since it is already divisible by 4. Thus $t\equiv3\bmod16$. The condition $4t\equiv1\bmod9$ now gives $t\equiv115\bmod144$. In bases 2 and 3 the dimension bounds at block length 4,096 are $0.5737$ and $0.5453$. Their sum is at most $1.119$. Taking $\epsilon=0.001$ in the two-base covering theorem gives $O(H^{0.12})$ uniformly.

In this last shape the three $y_i$ are $1,1,2$. Every prime $r\mid V$ is assigned to $C$ and is at least 5 since $V$ is odd and coprime to 3. Carry freedom of $V,V,2V$ in base $r$ implies that every digit of $V$ is at most $\lfloor(r-1)/4\rfloor$. Indeed the two copies of $V$ first add without a carry, yielding the digits of $2V$, and these double again without a carry. For $5\le r\le100000$ the dimension of the resulting digit set is at most $22/25$. This follows from $1+\lfloor(r-1)/4\rfloor\le(r+3)/4$, the increase of $(x+3)/x^{22/25}$ for $x\ge22$, and $100003^{25}<4^{25}100000^{22}$. The primes below 22 are checked individually.

For each such $r$ the three dimensions sum to at most $1.119+0.88<2$. Apply the moving-interval argument in the bases 2, 3 and $r$, using $t=(4V-1)/9$ at the third base. The finite union over these primes still gives $O_\epsilon(H^\epsilon)$. [Corso and Shmerkin, Corollary 1.17](https://arxiv.org/pdf/2409.04608)

Run `python3 check_paired_power_orders.py` with python-flint to reproduce the table, block counts and rational bounds, and `lean CarryOrder.lean` with Lean 4.12.0 to check the pairwise disjointness.
