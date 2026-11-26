# A recurring interval with two unrestricted prime-power factors

Let $J=10^{100}$ and $L=16384$. Throughout the interval
$$
 2^J\le t<2^{J+L},
$$
**Wasserman's full three-row conjecture holds for all the following weights.**

- $N=1+dt$ for every $1\le d\le1000$ when $t$ is a prime power.
- $N=1+5t$ and $N=1+7t$ when $t=p^a q^b$, for any primes $p,q$ and positive $a,b$.

**Intervals of this same multiplicative width recur with bounded gaps between their starting binary orders.** More precisely, fix any positive integer $D$. The set of $n\ge0$ such that the full conclusion holds with $J$ replaced by $J+Dn$ has lower natural density at least $2^{-6217}$ and has bounded gaps. The density bound is uniform in $D$. The gap bound may depend on $D$.

## Conditions that must survive

For the first family, a possible failure requires one of 16 combinations of slopes and forced bases from the [primitive denominator enumeration](shifted_directions.md) and the [subsequent exclusions](denominator_residue_hulls.md). Each requires that $(a t,b t,c t+1)$ add without carries in every one of its forced bases. For the second family, four cases suffice. Each requires that $kt,(d-k)t+1$ add without carries in the three listed bases.

| $d$ | $k$ | Bases |
|---:|---:|:---|
| 5 | 3 | 2, 3, 5 |
| 7 | 3 | 2, 3, 7 |
| 7 | 5 | 2, 5, 7 |
| 7 | 6 | 2, 3, 7 |

At even factor-seven weights, $U=6$ gives the last cut. At odd weights its three multipliers are $3,5,6$. Coincident primes or overlap with $d$ in this two-factor family leave at most two distinct predecessor primes, already settled. [Factor-five and factor-seven reductions](carry_intersections.md#necessary-additions)

## An interval certificate at each binary order

Fix $J\le j<J+L$ and normalize by $x=St/2^j$ with $S=2^{512}$, so $S\le x\le2S$ suffices. For each required base $p$ choose an exponent $E$ whose normalized root $v=S p^E/2^j$ is at least $2S$. An outward interval $[v_-,v_+]$ encloses the root. The checker uses certified logarithm bounds to choose $E$ initially, then encloses $p^E$ by repeated squaring with outward-rounded integer mantissas. On advancing $j$ it halves the bounds and multiplies by $p$ if needed. It restarts this computation every 4096 orders. At each odd prime it allows an extra leading digit and expands the interval to include $v(1+u)$ for every $|u|\le2^{-440}$. At 2 the root is exact.

A length-$k$ leading word of value $z$, including leading zeros to pad to $E$ digits, requires
$$
  v_- z/p^k\le x\le v_+(z+1)/p^k.
$$
The [carry graph](carry_intersections.md) tests the leading digits backward from the accepting states. Keep a stack of intersections of the indicated intervals with a state mask in each base. Subdivide a current word into its $p$ possible next digits, discarding only empty intersections and words with no predecessor states. Choosing the widest base interval affects only search order. With the interval deliberately closed on both sides, exhaustion proves that even the extra endpoints leave no admissible input. Reaching the digit or precision cutoff would return an unresolved result. All 327,680 cases exhausted their stacks, using 11,161,460 nodes in total.

A second checker uses Legendre's factorial floors directly. It obtains the leading powers from Arb ball arithmetic on $\exp(E\log p-j\log2)$ and descends through the powers in each base. Write the parts as $a_it+b_i$, their slope sum as $d$, and enclose a normalized layer $S p^h/2^j$ by $[l,u]$. (This checker uses $S=2^{640}$.) At a current left endpoint $x_0$, bound the total floor below by $\lfloor dx_0/u\rfloor$. Bound the part floors above by $f_i=\lfloor(a_ix_0+e_i)/l\rfloor$ where $e_i\ge b_i S/2^j$. If the total lower bound exceeds $\sum f_i$, a carry is certain from $x_0$ up to every $x$ for which $a_i x+e_i<(f_i+1)l$ for all $i$. Advance to the furthest certified right endpoint. The check succeeds only after these closed intervals cover $[S,2S]$. With the same relative margin at the odd primes, this independent cover used 2,352,590 certified intervals, with at most 1,310 in a single case. Both programs include small controls with exact multinomial or binomial coefficients, and both check their normalized powers against ordinary integer powers where these fit.

## Returning to the same leading powers

The 14 odd primes used across the 20 cases are
$$3,5,7,11,13,17,19,23,29,31,37,41,43,71.$$
For each of them put $\theta_p=\log_p2$. If $h\ge0$ satisfies $\|h\theta_p\|<2^{-444}$ at every such prime, the interval beginning at $J+h$ is certified by the same trees. Indeed, use $E'=E+\operatorname{round}(h\theta_p)$ at index $j+h$. The new normalized root differs from its old value by the factor
$$p^{\operatorname{round}(h\theta_p)-h\theta_p}.$$
Since $\log71<5$, its distance from 1 is less than $2^{-440}$. The enlarged bounds cover it and include enough leading digits for the whole binary interval. The new $E'$ is at least $E$ and the search finished within the original depth. At 2 the normalized root is identical and the extra offset in an affine part is still exactly the offset handled by the carry graph.

For an elementary density bound, partition the 14-dimensional unit torus into $K=2^{444\cdot14}$ half-open cubes of side $2^{-444}$. Among the $H+1$ points $mD\theta$ for $0\le m\le H$, at least $((H+1)^2/K-(H+1))/2$ unordered pairs lie in a common cube, by Cauchy–Schwarz. Each positive difference $n$ occurs in at most $H+1$ pairs. Every same-cube pair gives $\|nD\theta_p\|<2^{-444}$ at all primes. Thus at least $(H+1)/(2K)-1/2$ values in $[1,H]$ give certified starts, which proves the density claim without an independence hypothesis on the logarithms.

For bounded gaps, the closure of the nonnegative orbit of $D\theta$ is a compact subgroup: pigeonholing gives positive multiples approaching zero and hence negative multiples in the closure as well (the finite-orbit case is immediate). Every forward orbit in this subgroup is dense there. If $U$ is the open neighborhood given by the 14 strict inequalities, the translates $U-kD\theta$, $k\ge0$, cover the subgroup. A finite subcover shows there is a visit within a bounded number of steps after every start.

## Reproduction

Run from the repository root with Python and python-flint. The prefix checker uses standard-library arithmetic for every inequality after the enumeration of the possible rows. The second checker independently computes the powers and the exclusion intervals.

~~~~sh
python3 check_distant_orders.py --include-cuts
python3 check_distant_orders_cross.py --include-cuts
~~~~
