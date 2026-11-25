# Odd sparse weights and the endpoint

**Theorem.** Fix $s\ge3$. Among odd weights $N\le X$ having exactly $s$ binary ones, at most $O_s((\log X)^{s-2})$ can admit a counterexample. The total number of odd weights in this class has order $(\log X)^{s-1}$. More precisely, outside a finite set, any exceptional exponent vector satisfies one of finitely many proper affine linear equations on the differences between its positive bit positions.

The [bounded-digit theorem](bounded_binary_weight.md) reduces large weights in this family to the endpoint $E=(1,1,N-2)$. For $s\le5$ it already excludes a counterexample at every sufficiently large odd weight. The following argument for $s\ge6$ bounds the remaining exceptions without imposing an upper limit on $s$.

Write $N=1+Qm$, where $Q=2^k$, $k\ge1$, and $m$ is odd. A row avoiding 2 has the form
$$A=(Qa,Qb,Qc+1),\qquad m=a+b+c,$$
where the bits of $a,b,c$ partition those of $m$. In a selection containing $E$, the two other rows are imprimitive, so $a,b,c>0$. Assign the prime powers of $N(N-1)$ to these two rows as in the [endpoint reduction](endpoint_reduction.md), assigning 2 to $A$, which is possible by construction. Write $U,V$ for their products from $N$ and $S,T$ for their products from $m$, respectively. Thus $UV=N$ and $ST=m$. The second row has parts $Vy_i$ with $\sum y_i=U$ and $T^2\mid y_1y_2y_3$. We obtain
$$U\mid\gcd(a,b),\qquad S\mid\gcd(ab,m),\qquad 27(m/S)^2\le U^3. \tag{1}$$
For the middle claim, any full power assigned to $S$ divides two of the three actual parts of $A$, hence divides at least one of $Qa,Qb$. It is odd.

There are $s-2$ variable bit positions in $m$ besides position zero. Treat their powers of two as the coordinates of a torus. For each fixed bit partition, $a,b,c$ are sums of disjoint sets of these variables and possibly the constant 1. The polynomial pairs $(a,b),(a,m),(b,m)$ are coprime: for example a factor common to $a,m$ would divide the nonzero polynomial $b+c$ in disjoint variables. Constant cases give the same gcd bound directly. Use the outside-$S$ gcd theorem with places $2,\infty$ and exponent $1/10$ on the three pairs. The coordinate heights are at most $\log m$. Off finitely many translated proper subgroups, for all sufficiently large $m$ the odd parts of all three gcds are at most $m^{1/10}$. Thus $U\le m^{1/10}$ and $S\le m^{2/10}$. Equation (1) would give $m=ST\le S(U/3)^{3/2}< m^{7/20}$, impossible. The prime powers assigned to $U$ are all odd here. [Levin–Wang, Theorem 1.13(a)](https://doi.org/10.1515/crelle-2019-0033)

There are finitely many partitions. Each relevant translate imposes a nontrivial affine linear equation on the positive bit positions of $m$, as in the bounded-digit proof. Bounded values of $m$ impose fixed positions too. The bit positions of $m$ are the positive bit positions of $N$ shifted by $-k$. Lifting the equations by this substitution gives proper affine hyperplanes in the $s-1$ positive exponent coordinates of $N$. A nonzero coefficient allows one coordinate to be recovered from the rest, giving the asserted count. In fact (1) excludes the endpoint for **every** $k\ge1$ at a multiplier outside the exceptional set.

## A necessary pairing when $s=6$

At odd six-bit weights the only possible binary row accompanying $E$ divides the six bits into three pairs. A singleton would force that row to be primitive. Let $p$ be an odd prime factor of the common divisor of the three parts. Each difference $j$ between the two bit positions of a pair satisfies $2^j\equiv-1\pmod p$. Consequently all three differences have the same 2-adic valuation $t=v_2(\operatorname{ord}_p(2))-1$. Each pair matches an exponent in one class modulo $2^{t+1}$ to an exponent in the class displaced by $2^t$. Equivalently, the bit-position polynomial $\sum_{2^i\text{ occurring in }N} X^i$ must be divisible by $X^{2^t}+1$. This is a necessary condition at every weight. For example, if all five positive bit positions are odd, then the pair containing position zero forces $t=0$, but the other two pairs cannot match even to odd. All sufficiently large such weights therefore satisfy the full conjecture.

## Every odd six-bit weight below $2^{64}$

The full conjecture holds at all $\binom{63}{5}=7{,}028{,}847$ such weights. Write $N=1+Qm$ as above. In an endpoint candidate each part of $A$ uses two bits, so $a,b$ use two bits of $m$ apiece and $c$ is its remaining single power of two. This gives 15 assignments, including the choice of which part contains the final 1.

Equation (1) first tests each assignment independently of $k$, using the odd part of $\gcd(a,b)$ and $\gcd(ab,m)$ as upper bounds divisible by $U,S$. For a particular $k$, both $U$ and $S$ also divide $Qc+1$. This is immediate for $U$. For a full power assigned to $S$, if both $a,b$ were divisible by its prime then so would be $c=m-a-b$, impossible for a power of two and an odd prime. Thus the two vanishing positions must include the shifted part.

There were 8,635 binary rows above the one-million cutoff requiring further checks. Factor $N$ and $m$ and use factorial valuations on $A$. For each full power from $N$ or $m$ at whose prime $A$ has positive valuation, the last row has to avoid that prime. The [exact CRT row search](endpoint_reduction.md) checks these remaining conditions, including the higher base-$p$ digits. It finds no row satisfying them. The input factorizations are verified by prime proofs and multiplication.

Of the 7,028,847 weights, 11,627 are covered by the [one-million theorem](full_range_1m.md). For the case without $E$, 6,968,699 more are excluded by the gap and height-six bound from the [four- and five-bit verifier](bounded_binary_weight.md). At the remaining 48,521 weights, test all 90 partitions of six bits into three nonempty sets. The small-sum test excludes 24,178 rows, the two-ordinary-row product excludes 4,196,472, and the near-prime product excludes 146,240. These counts sum to $90\cdot48{,}521$.

```sh
python3 check_six_binary_digits.py 64
```

## An unbounded-exponent check for 595,665 multipliers

**Theorem.** Let $m<2^{64}$ be odd with exactly five binary ones. For **every** $k\ge1$, the endpoint is impossible at $N=1+m2^k$. Consequently the full conjecture holds at every sufficiently large weight in this union of 595,665 families. The last implication uses the [bounded-digit theorem](bounded_binary_weight.md); the endpoint conclusion has no threshold.

Write $Q=2^k$ and use the 15 normalized binary rows above. For each assignment from (1), enumerate every $U\ge3$ dividing the odd part of $\gcd(a,b)$ and coprime to $m$. Choose $S$ as any product of full powers of $m$ each dividing $a$ or $b$, and put $T=m/S$. The remaining row is $(Vy_1,Vy_2,Vy_3)$ with $UV=N$ and $\sum y_i=U$. For each full power $w\mid T$ assigned to it, two of its normalized parts are divisible by $w$, and the third is $U\pmod w$. Thus $2w+1\le U$ and $27T^2\le U^3$. Assign the exceptional positions among the three bins, take least positive CRT representatives modulo $T$, and add all nonnegative lifts of sum $U$. This enumeration makes no choice of $k$.

Of 29,093 preliminary assignments passing the initial product inequality, the scale and CRT enumeration yields just 726 normalized row pairs. Two one-digit checks eliminate 725 of them.

- For a full power $w=p^j$ of $U$, the $j$-th base-$p$ digits of the two unshifted parts are $(a/w)Q\bmod p$ and $(b/w)Q\bmod p$. Because $U\mid Qc+1$, here $Q\equiv-c^{-1}\pmod p$. Their sum forces a carry for 677 pairs.
- For a full power $w=p^j$ of $T$, the $j$-th digits of the two parts divisible by $w$ are $(y_i/w)V\bmod p$. We know $V\equiv U^{-1}\pmod p$, as $UV=1+mQ$. This forces a carry for 48 further pairs.

For the final pair,
$$m=6917529027643179011,\quad U=768614336404564651,\quad S=1,\quad c=2^{21},$$
and the sorted $(y_i)$ are $(167885246399975556,238535483711833759,362193606292755336)$. The condition $U\mid 2^k c+1$ gives $k\equiv40\pmod{122}$. Take $p=23\mid T$. As $k$ traverses this class there are 11 residues for $Q$ modulo 23. Each determines $V=(mQ+1)/U$ modulo $23^2$, since $23\mid m$. For every one, the sum of the three residues of $Vy_i$ modulo $23^2$ is at least $23^2$, forcing a carry. This finishes the exhaustion.

The standard-library integer and CRT steps use the prime-proved factorizations from python-flint. Run the verifier from the repository root.

```sh
python3 check_six_bit_rays.py 64
```
