# Shifted proportional rows at binary prefixes

Use the [determinant theorem](projected_rows.md) with $N=d2^k+c$, $1\le c<2^k$ and $(dc)^5<27(N-1)^2$. One may absorb any power of two in $d$ into $2^k$ before testing these conditions. Here take $d$ odd.

**Theorem.** For odd $d\le4095$ other than $1887,3451,4095$, every counterexample requires the *unshifted* proportionality $b=mz$ in the binary row. In particular the conjecture holds whenever $\gcd(d,c)\le6$ at such a prefix. For the three specified prefixes this conclusion holds through $N\le10^{1000}$. Over unbounded weights at most $O_\epsilon(X^\epsilon)$ distinct $N\le X$ can use the shifted case at these prefixes, for each $\epsilon>0$.

Unshifted proportionality requires $s\mid\gcd(d,c)$ with an odd $s\ge7$. The theorem removes the separate test on $\gcd(d,c-1)$ at the indicated prefixes, even when $d\mid c-1$.

## Counting exponents at an unrestricted prefix

For every fixed $d$ there is a $\theta_d<1$ such that, for any fixed positive offset $c$ and any interval of $K\ge2$ consecutive positive exponents, at most $O_d(K^{\theta_d})$ exponents can meet the hypotheses and use a shifted row. The constant is uniform in $c$ and the location of the interval. For odd $d\le4095$ one may take $\theta_d=0.34$.

Indeed, in the shifted case the first row is
$$A=(xt,yt,zt+1),\qquad t=g2^k+m,\qquad s=x+y+z,\qquad N=st+1.$$
Here the upper parts are $g(x,y,z)$, the integers $x,y,z$ are positive and primitive together, and $gs=d$. The second row has the form $V(f_1,f_2,f_3)$ with $UV=N$ and $\sum f_i=U$. The assigned scale $U\ge3$ divides $x,y$ and is coprime to $t$ and $s$. The choices of $g,x,y,z,U,f_i$ are all bounded in terms of $d$, and $c=sm+1$ fixes $m$ for a given choice.

If $U$ has an odd prime factor $r$, the first row must avoid $r$. Otherwise $U$ is a power of two, and choose an odd prime $r\mid s$, possible since $s$ is odd and at least 7. At least one of the two rows must avoid $r$. In the first row two positive multiples $xt,yt$ must add without a carry in that base. In the second row the same holds for $f_1V,f_2V$.

For any positive multipliers at most $d$, a run of $L=1+\lceil\log_r d\rceil$ consecutive digits $r-1$ in the input forces the last digit of both products on that block to be $r-1$. The incoming multiplication carry for a factor $a$ lies between 0 and $a-1$. The first $L-1$ maximal digits take it to $a-1$, as is seen by dividing $a(r^{L-1}-1)$ plus the incoming carry by $r^{L-1}$. Thus every admissible input omits this block.

As $k$ varies in the given interval, any residue of $t=g2^k+m$ or $V=(sg2^k+sm+1)/U$ modulo $r^l$ occurs $O_d(1)$ times when $r^l\asymp_d K$. For the second expression $U$ is a power of two, invertible modulo $r^l$. After cancelling the fixed $r$-part of the slope this follows from $\operatorname{ord}_{r^n}(2)\gg_r r^n$. In detail, if $h$ is the order modulo $r$ and $v=v_r(2^h-1)$, then $v_r(2^{hb}-1)=v+v_r(b)$ by binomial expansion at the odd prime $r$. Choose $l$ a multiple of $L$ and count allowed residues in $l/L$ disjoint blocks. The exponent $\log(r^L-1)/(L\log r)$ is below one. Maximizing over the finitely many possible bases proves the assertion. This uses the same order calculation as the [prime-power ray estimate](prime_power_ray_density.md).

## Finite classification of the row shapes

Enumerate positive upper parts $a,b,c'$ with pairwise disjoint bits and odd sum at most 4095, taking $a<b$ for the two unshifted positions. Divide by their common gcd to get $x,y,z$. For each divisor $U\ge3$ of $\gcd(x,y)$ use $t\equiv-s^{-1}\pmod U$ and $\gcd(U,V)=1$.

At a full power $w=r^e$ of $s$, if $r\mid t$ then $A$ cannot avoid $r$. Its shifted part is 1 modulo $r$, so both $x$ and $y$ would have to be divisible by $w$ to reach the full power in $st$. This would also make $r\mid z$, contradicting primitivity. If $r\nmid t$, avoidance by $A$ requires one of $x,y$ to be divisible by $w$ and the other to be a unit at $r$. Both cannot be divisible by $r$ since $r\mid s$. In this case precisely one of the two unshifted parts is nonzero modulo $w$, so $(x+y)t\equiv1\pmod w$.

Every power $w$ failing this test must go to the second row. Two of its normalized parts must be multiples of $w$, so $2w+1\le U$. If the product of required powers is $T_0$, then $27T_0^2\le U^3$. These size tests leave 3,976 distinct choices of $(s,x,y,z,U)$. Test the first row in base two and the bases of $U$. Add any further powers of $s$ whose corresponding carry test on $A$ is impossible to the list required in the second row. Enumerate all second rows modulo the required powers, including positive lifts, and test their carry conditions. Whenever two normalized parts cannot be zero at a power of $s$, impose its congruence on $t$ and test the resulting combined progression in all bases forced on $A$.

For an affine progression one can test a whole prime base without bounding $t$. Read its parameter from the low digits upward and track the carry of each affine part. A digit is allowed only if the output digits sum without a carry. Reachability of a state whose terminal carries add freely decides whether any nonnegative parameter works. At a prime dividing $U$, the very first digit also imposes $r\nmid V$. The [carry-graph construction](carry_intersections.md) gives the transitions. Tests at individual bases can only remove a shape when an entire required language is empty.

Of the 3,976 choices, 360 pass the initial tests on $A$, and 254 pass the refined size tests. These produce 23 unordered complementary assignments before the remaining carry tests. Generating all admissible single parts and testing their sums gives exactly the same rows as the CRT assignment in every one of the 254 cases. Through prefix 511 an independent enumeration starting with primitive directions also reproduces all eight preliminary choices. The three survivors are as follows. Each primitive slope has pairwise disjoint bits. Its positive odd multiple cannot give another prefix through 4095 since the next odd multiplier is at least three.

| $s=d$ | $(x,y,z)$ | $U$ | $(f_1,f_2,f_3)$ | Necessary progression for $t$ | Bases forced on $A$ |
|---:|---|---:|---|---|---|
| 1887 | (833, 1036, 18) | 7 | (1, 3, 3) | $1223\pmod{4403}$ | 2, 7, 17, 37 |
| 3451 | (1281, 2074, 96) | 61 | (3, 29, 29) | $6276\pmod{7259}$ | 2, 7, 17, 61 |
| 4095 | (1946, 2085, 64) | 139 | (9, 13, 117) | $76\pmod{4865}$ | 2, 5, 7, 139 |

For each row, intersect the required carry languages with its progression through $t\le10^{1000}$. All three intersections are empty. The interval traversal takes respectively 144,968, 153,826 and 160,003 nodes. The tests start at $t=2$, and $t=g2^k+m\ge2$.

A separate verifier covers the same progressions by intervals obtained directly from factorial floors. At a progression point $t$, let the parts be $a_i= z_i t+\delta_i$. If $\lfloor N/w\rfloor>\sum_i\lfloor a_i/w\rfloor$ at a power of an assigned prime, the valuation is positive throughout $[t,t+B]$ for
$$B=\min_i\left\lfloor\frac{(\lfloor a_i/w\rfloor+1)w-1-a_i}{z_i}\right\rfloor.$$
All three denominator floors stay fixed while the numerator floor can only increase. Advance to the first progression point beyond the largest such certified endpoint. This verifies the three entire ranges in 8,471, 8,736 and 5,973 steps.

For an unbounded count use 64-digit blocks in two of the forced bases. Count paths from every reachable starting state, with any terminal state. The corresponding integer missing-digit set bounds are 32% in base 2 and 34% in base 7 at $s=1887$, 34% in each of bases 2 and 7 at $s=3451$, and 32% in base 5 and 33% in base 2 at $s=4095$. Each is checked by the exact inequality $T^{100}\le r^{64j}$ for the path count $T$ and the listed percentage $j$. The base-$r^{64}$ digits allowed by the paths include every aligned block of an accepted $t$, including leading zero padding. The two dimensions sum to less than one. Theorem B of Glasscock, Moreira and Richter therefore gives $O_\epsilon(X^\epsilon)$ possible parameters through $X$, as in the [earlier intersection argument](carry_intersections.md). The bound on just the odd-base blocks, combined with the multiplicative-order argument above, gives $O(K^{0.34})$ exponents in an interval of length $K$. [Glasscock, Moreira and Richter](https://arxiv.org/abs/2007.05480)

Run `python3 check_shifted_prefixes.py` for the classification and both independent interval checks. The default uses 12 prefix bits and the bound $t\le10^{1000}$.
