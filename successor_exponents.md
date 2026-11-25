# Successor weights with the two exponents in either order

Let $N=p^a(p^j+d)+1$, where $p$ is prime and $a,j,d$ are positive integers. Wasserman's three-row conjecture holds for **all $a,j$ and every $1\le d\le191$**. More precisely, through $d=256$ it holds with the same unrestricted exponents except possibly on the five rays below, and it holds on those rays as well for **every $a\le10^{12}$**.

| $d$ | $p$ | $j$ | Necessary exponent class |
|---:|---:|---:|:---|
| 192 | 1459 | 1 | $a\equiv86\pmod{126}$ |
| 198 | 1453 | 1 | $a\equiv20\pmod{126}$ |
| 224 | 1427 | 1 | $a\equiv2\pmod{126}$ |
| 237 | 4889 | 1 | $a\equiv391\pmod{870}$ |
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

Through $D=256$ there are 14,007 row pairs with height and prime choices after the power-orbit test. For each base take the largest power at most 256, or the base itself if larger. The simultaneous checks retain 127 cases. Repeating with cap 4,096 leaves 56. For these 56, put $P=R+Mn$ and test whether any nonnegative integer $n$ makes each forced affine row carry-free in its required base. The [finite carry automaton](carry_intersections.md) reads the digits of $n$ and exhausts all reachable multiplication-carry states. Different bases can use different $n$ in this necessary test. It leaves five cases, precisely those in the table.

At $Q=1651$ the only surviving shapes are $U=4$, $T=1$, $(\alpha,\beta,\gamma)=(104,1524,23)$ and $(f,g,h)=(1,1,2)$, with $P\equiv2297\pmod{6604}$. At $Q=5126$ they are $U=7$, $T=2$, $(231,4893,2)$ and $(1,2,4)$, with $P\equiv8970\pmod{17941}$. Computing the orders and discrete logarithms by a complete orbit gives the exponent classes in the table.

## The remaining five exponent intervals

Each case forces an addition in base two. Begin with the necessary exponent class $a_0$ modulo $L$ and lift it at the successive layers $w=2,4,8,\ldots$. If $p^L\equiv1\pmod w$, every representative in the class still has the same power residue. Otherwise split it into the two classes modulo $2L$ and double the stride. This suffices because the old stride gave 1 at the preceding layer and its square gives 1 at the new layer. Discard classes failing the required binary addition at $w$, and those with least representative exceeding $10^{12}$. Once the stride exceeds $10^{12}$ each remaining class represents at most one exponent in the interval. For every such exponent some required base has a carry in its first 100 digits. All powers in this computation are evaluated by modular exponentiation.

Run the following with Python 3 and python-flint.

```sh
python3 check_successor_exponents.py
python3 check_successor_exponents_cross.py
```

The cross-check independently constructs the power residues modulo $U B$ and enumerates both normalized rows directly at three selected assignments. It agrees on all five row pairs there. It also compares 343 generalized CRT examples with exhaustive search over their full periods, compares binary lifting with direct powers, and includes a positive base-two control checked by the coefficient $\binom{244}{36}\binom{208}{144}$.
