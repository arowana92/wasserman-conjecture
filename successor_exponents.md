# Successor weights with the two exponents in either order

Let $N=p^a(p^j+d)+1$, where $p$ is prime and $a,j,d$ are positive integers. Wasserman's three-row conjecture holds for **all $a,j$ and every $1\le d\le191$**. More precisely, through $d=256$ it holds with the same unrestricted exponents except possibly on the four rays below, and it holds on those rays as well for **every $a\le10^{24}$**.

| $d$ | $p$ | $j$ | Necessary exponent class |
|---:|---:|---:|:---|
| 192 | 1459 | 1 | $a\equiv86\pmod{126}$ |
| 198 | 1453 | 1 | $a\equiv20\pmod{126}$ |
| 224 | 1427 | 1 | $a\equiv2\pmod{126}$ |
| 249 | 4877 | 1 | $a\equiv587\pmod{870}$ |

These are possible locations after sieving. The proof uses necessary conditions for a pair of rows to avoid all primes in the endpoint coefficient. For the ordered exponents $j\ge a$, the [paired theorem](successor_powers.md) already treats $d\le1000$.

## Enumerating normalized rows without bounding the exponent

Use $P=p^a$, $H=p^j$, $Q=H+d$, and the notation $U,B,T,x,y,t,e$ in the [proof of the size bound](successor_powers.md). For any $a,j$ it gives $27Q^2<d^5$, $Q=BT$, $27T^2\le U^3$, $3\le U\mid y<d$, $e\le d$ and $B/e\mid y/U$. Also $e$ is coprime to $y$, since $e\mid Py-1$. This determines $e=B/\gcd(B,y/U)$. Every full power in $Q$ is at most $d$, by the bounds on $e$, $B/e$ and the individual powers of $T$.

For a maximum offset $D$, enumerate $Q\le\lfloor\sqrt{D^5/27}\rfloor$ and all prime powers $H$ in $[Q-D,Q-4]$. Remove the full $p$-power from $Q$, obtaining $q_0$. The [cofactor theorem](prime_power_cofactors.md) and the [all-weight endpoint reduction](all_weights_endpoint_reduction.md) already settle $q_0\le1000$ except for the fifteen listed cofactors with $p$ odd. For the remaining $Q,H$, enumerate every unitary $T\mid Q$ coprime to $p$, every admissible $U$, and $y$. The other small part satisfies either $x\equiv0\pmod{Ue}$ or $x+y\equiv0\pmod e$ and $x+y\equiv Q\pmod U$. Retain $x+y\le Q-H$. The corresponding first row has shape
$$A=(\alpha P,\beta P,\gamma P+1)$$
with $(\alpha,\beta,\gamma)$ adding to $Q$ without carries in base $p$.

We must have $P\equiv-Q^{-1}\pmod U$. For a full power $w$ in $B$, two unshifted slopes divisible by $w$ need no additional congruence. Otherwise one must vanish modulo $w$ and the other must be a unit, and $P(\alpha+\beta)\equiv1\pmod w$. Indeed the three actual residues are $0,0,1$. At the base $p$ the full power in $PQ$ forces the normalized slopes to be zero modulo the $p$-part of $Q$. Combining the congruences gives $P\equiv R\pmod M$ with $U\mid M$ and $R$ a unit. Keep only $R$ in the multiplicative orbit of $p$ modulo $M$. Independently enumerate the complementary row $C=V(f,g,h)$, where $V=(QP+1)/U$ and $f+g+h=U$. For each full power $w$ in $T$, the three slopes modulo $w$ are $0,0,U$ up to permutation. Enumerating by exceptional positions and by admissible individual parts gives identical sets for every call.

## Combining prime-power orbits across bases

The first row must avoid the bases of $UB$, and the second must avoid those of $T$. If $p,Q,U$ are all odd, $2\mid V$ so the second row must also avoid 2. The first row at $p$ has already been tested through its normalized digits.

Fix one of these required bases $r\ne p$. At modulus $w=r^k$, all the low digits of the required row follow from $v=p^a\bmod w$. For the second row use $V\equiv(Qv+1)U^{-1}\pmod w$. Require the three residues to sum below each layer $r^i$, $1\le i\le k$. For $r\mid U$, if $w$ reaches the layer $r^{v_r(U)+1}$ then $Qv+1$ must be nonzero at that layer, by unitarity. Each such check gives a set of allowed exponents modulo the order of $p$ at $w$. Include the class from $p^a\equiv R\pmod M$ as another condition.

To test simultaneous compatibility for periods $m_i$, project the allowed residues modulo $m_i$ to $\gcd(m_i,\operatorname{lcm}_{j\ne i}m_j)$. These projections retain every pairwise gcd modulus. The generalized Chinese remainder theorem therefore makes compatibility of representatives equivalent before and after projection. The program searches the projected sets by exact integer congruences.

Through $D=256$ there are 14,007 row pairs with height and prime choices after the power-orbit test. For each base take the largest power at most 256, or the base itself if larger. The simultaneous checks retain 127 cases. Repeating with cap 4,096 leaves 56. For these 56, put $P=R+Mn$ and test whether any nonnegative integer $n$ makes each forced affine row carry-free in its required base. The [finite carry automaton](carry_intersections.md) reads the digits of $n$ and exhausts all reachable multiplication-carry states. Different bases can use different $n$ in this necessary test. It leaves five cases before the fixed-divisor test below.

At $Q=1651$ the only surviving shapes are $U=4$, $T=1$, $(\alpha,\beta,\gamma)=(104,1524,23)$ and $(f,g,h)=(1,1,2)$, with $P\equiv2297\pmod{6604}$. At $Q=5126$ they are $U=7$, $T=2$, $(231,4893,2)$ and $(1,2,4)$, with $P\equiv8970\pmod{17941}$. Computing the orders and discrete logarithms by a complete orbit gives the exponent classes in the table, together with $a\equiv391\pmod{870}$ for $p=4889$.

## A common divisor on the fifth ray

For $p=4889$ the last congruence forces $a$ odd. Since $Q=5126\equiv1\pmod5$ and $4889\equiv-1\pmod5$, the scale $V=(Qp^a+1)/7$ is divisible by 5. Thus the complementary row must avoid 5. But in $V(1,2,4)$ its coefficient is always divisible by 5. At the first nonzero base-five digit of $V$, the residues have sum $7,9,6,8$ when the digit is $1,2,3,4$, respectively. This eliminates the ray for every exponent and settles offset 237 in full.

## The four remaining exponent intervals

At $Q=1651$ the first row must avoid 2, and at $Q=5126$ it must avoid 7. Use that base $r$ to lift the exponent class $a_0$ modulo $L$ through $w=r,r^2,r^3,\ldots$. Here $p^L\equiv1\pmod r$ at the outset. If $p^L\equiv1\pmod w$, the class still has a single power residue. Otherwise split it into the $r$ classes modulo $rL$ and multiply the stride by $r$. This suffices because the old stride gave 1 at the preceding layer and its $r$-th power gives 1 at the new layer. Discard classes failing the required addition at $w$, and those with least representative exceeding $10^{24}$. Once the stride exceeds $10^{24}$ each remaining class represents at most one exponent in the interval. For every such exponent some required base has a carry in its first 100 digits. All powers in this computation are evaluated modulo integers of bounded size.

Run the following with Python 3 and python-flint.

```sh
python3 check_successor_exponents.py
python3 check_successor_exponents_cross.py
```

The cross-check independently constructs the power residues modulo $U B$ and enumerates both normalized rows directly at three selected assignments. It agrees on all five row pairs there. It also compares 343 generalized CRT examples with exhaustive search over their full periods, compares all surviving base-two and base-seven lifts with direct powers through exponent one million, and includes a positive base-two control checked by the coefficient $\binom{244}{36}\binom{208}{144}$.

## Counting possible exceptions at any fixed offset bound

Fix $D\ge1$. Among **all** weights $p^a(p^j+d)+1\le X$ with $1\le d\le D$ and **odd** $p$, the number of possible exceptions is
$$O_D\!\left((\log X)^{\log_2\varphi}\right),\qquad \varphi=(1+\sqrt5)/2.$$
Including $p=2$, it is $O_D((\log X)^{1-\eta_D})$ for some $\eta_D>0$.

Indeed (1) in the [size-bound proof](successor_powers.md) leaves finitely many $(p,H,Q)$ depending only on $D$, including for $a>j$. The bounds $U<d$ and $T\mid Q$ also leave finitely many normalized row shapes. For odd $p$, at least one of $A,C$ must avoid 2. When $A$ does, the binary expansions of $\alpha P$ and $\beta P$ are disjoint. When $C$ does, two specified positive multiples of $V$ are disjoint. The [binary digit encoding](projected_rows.md#counting-entire-families-of-offsets) bounds the number of low $k$-bit scalars obeying either condition by $O_D(\varphi^k)$. The encoding reads digits from low to high, so this bound applies to their residues as well. Cancelling common powers of two loses only a bounded number of top digits. For equal valuations the two odd multiples have the same lowest one, so disjointness at the low layers forces all remaining input digits there to be zero.

For the $C$ condition it suffices to consider $U$ odd, since $2\mid U$ already forces $A$ to avoid 2. The map $P\mapsto (QP+1)U^{-1}\pmod{2^k}$ then has at most $2^{v_2(Q)}$ preimages per value. For fixed odd $p$, the order of $p$ modulo $2^k$ is bounded below by a positive constant times $2^k$, by the [order-lifting argument](prime_power_ray_density.md#counting-exponents). Taking $2^k$ comparable to $\log X$, each allowed residue represents only $O_D(1)$ exponents in the interval. This proves the first count. For $p=2$, absorb the two-part of $Q$ into the exponent and use the [power-saving cofactor-ray bound](prime_power_ray_density.md) for each of the finitely many cofactors, together with the endpoint reduction.

At $D=256$ the four explicit rays give a sharper bound. At $Q=1651$, the binary carry graph for $(104P,1524P,23P+1)$ has 152 states. The total number $Z$ of readable 1,024-digit paths from all states satisfies $Z^5<2^{1024}$. At $Q=5126$ the base-seven graph for $(231P,4893P,2P+1)$ has 742 states and $Z^{16}<7^{5\cdot1024}$. Both comparisons use integer recurrence in the cross-check script. Splitting low-digit strings into 1,024-digit blocks gives $O(2^{k/5})$ binary strings or $O(7^{5k/16})$ base-seven strings of length $k$. In each case the order of the fixed prime $p$ modulo $r^k$ is bounded below by a constant times $r^k$. Thus the exceptional weights through $X$ number $O((\log X)^{1/5})$ on the three $Q=1651$ rays and $O((\log X)^{5/16})$ on the single $Q=5126$ ray.
