# A single perfect power when one predecessor power divides the other minus one

The [leading-digit theorem](nested_predecessor_cofactors.md) now settles the full nested family at cofactor ten for all exponents.

Let $P<Q$ be powers of distinct primes with $P\mid Q-1$. At $N=1+10PQ$, the full trinomial conjecture holds **if $P$ is a square or if $Q$ is a fourth power**. There is no bound on the other prime power or on the digits of either. This complements the [sparse and even cases](ten_nested_powers.md).

More generally let $b$ be a prime power. At $N=1+bPQ$ it suffices that either $P$ or $Q$ is one modulo $b$. At $N=1+2bPQ$ with $b$ odd, it suffices that either varying power is itself a perfect $e$-th power for any positive multiple $e$ of
$$
 \operatorname{lcm}\big(\varphi(b),\ \ell-1:\ \ell<b\text{ an odd prime}\big).
$$
At $1+bPQ$ the corresponding period is just $\varphi(b)$. Either statement permits shared prime factors with the static cofactor.

## Two congruences at each of two moduli

Write $D=b$ or $2b$ as above. The [static scale theorem](static_square_cofactors.md) and [distinct exceptional positions](predecessor_inverse_digit_gaps.md) show that any failure must have $(b,PQ)=1$ and a row taking the prime bases of $b,P,Q$. For $Q=mP+1$ the [row calculation](ten_nested_powers.md#determining-the-candidate-row) works for general $D$ and gives
$$
 \big(xPQ,\ (tP-1)Q+1,\ (kP+1)Q\big),\qquad x+t+k=D,\quad t,k\ge1.
$$
In detail, the second part is a multiple of $P$ congruent to one modulo $Q$, and all such parts are $P(tQ-m)=(tP-1)Q+1$. The third part cannot equal $Q$, since the row scale $U\ge3$ is coprime to $Q$. We have $U\mid x$, $U\mid N$, and both of the last two parts are divisible by $b$. Therefore
$$
 kP\equiv-1\pmod{bU},\qquad xQ\equiv-k\pmod b,\qquad DQ\equiv k\pmod U. \tag{1}
$$
Here $(U,DbPQ)=1$ and $(x,b)=1$. The second congruence follows from $xPQ\equiv1\bmod b$; the third follows by substituting $P^{-1}\equiv-k\bmod U$ into $DPQ\equiv-1\bmod U$.

For $D=b$, $P\equiv1\bmod b$ would force $k\ge b-1$, and $Q\equiv1\bmod b$ would force $x+k\ge b$. Both contradict $x+k<D$. More generally, write $[z]_b$ for the least nonnegative residue. Necessarily
$$
 [(PQ)^{-1}]_b<[P^{-1}]_b,
$$
since $x=[(PQ)^{-1}]_b$ and $k=b-[P^{-1}]_b$. For each fixed unit $Q\not\equiv1\bmod b$, exactly half the unit classes of $P$ fail this strict inequality: negating $P$ reverses it, and equality would force $Q\equiv1\bmod b$.

One consequence is a criterion symmetric in three prime powers. If $A<B<C$ are powers of distinct primes, the conjecture holds at $1+ABC$ whenever **at least two** of $A\mid B-1,\ A\mid C-1,\ B\mid C-1$ hold. If $B\mid C-1$, use the preceding result with static cofactor $A$. In the remaining case $BC\equiv1\bmod A$ and the [static numerator theorem](static_square_cofactors.md#at-most-two-possible-numerators) would require its numerator to be 1 and at least 3.

For $D=2b$, first suppose $P\equiv1\bmod b$. Then $k=b-1$ and $x\le b$. As $x$ is a unit modulo $b$, $U\le x<b$. Choose any prime $\ell\mid U$. It is odd and coprime to $bP$. If $P$ is a perfect power with the indicated exponent, it is one modulo $\ell$ and modulo $b$. Now $kP+1\equiv b\not\equiv0\bmod\ell$, contradicting (1).

If instead $Q\equiv1\bmod b$, then $x+k=b$, so again $U<b$. The last congruence in (1) gives $b(2Q-1)\equiv0\bmod U$ since $U\mid x$. With $Q$ a perfect power of the stated exponent, it is one modulo $b$ and each prime of $U$, a contradiction. These proofs use the power condition on one side only.

## Smaller exponents

For either side, replace $U$ in (1) by any odd prime factor $u$, or by $u=4$ if $U$ is a power of two. Thus
$$
 u\mid x,\quad (u,Dk)=1,\quad x+k<D.
$$
If $P$ is an $e$-th power, $-1/k$ must be a unit $e$-th power modulo both $b$ and $u$. If $Q$ is an $e$-th power, the corresponding residues are $-k/x$ modulo $b$ and $k/D$ modulo $u$. The table gives every $b\le500$ excluded by these necessary tests. Entries are the **static prime powers $b$**, so the weight in the last two columns is $1+2bPQ$.

| Exponent | $D=b$, power $P$ | $D=b$, power $Q$ | $D=2b$, power $P$ | $D=2b$, power $Q$ |
|---:|---|---|---|---|
| 2 | 2,3,4,5,7,8,9,16 | 2,3,4,5,8 | 3,5 | 3 |
| 3 | 2,3,4 | 2,3,4,7,9 | 3 | 3 |
| 4 | 2,3,4,5,7,8,9,16,32 | 2,3,4,5,8,9,16 | 3,5 | 3,5 |
| 6 | 2,3,4,5,7,8,9,16 | 2,3,4,5,7,8,9,13,16 | 3,5,9 | 3,7,9 |

At ten there is also a short check. The nonempty choices of $(U,x)$ are $(3,3),(3,6),(7,7)$. For $U=3$ and square $P$, the units modulo 15 that are squares are 1 and 4. Thus $k\equiv14$ or $11\bmod15$, impossible for $1\le k\le6$. For $U=7$, $k=1,2$ and $kP\equiv-1\bmod7$ would make a square times a square into a nonsquare. For fourth-power $Q$, (1) at $U=3$ gives $k\equiv1\bmod3$. Modulo 5 it gives $k\equiv-x\bmod5$ since $Q\equiv1\bmod5$. At $x=3$ neither eligible $k=1,4$ works; at $x=6$ the sole eligible $k=1$ fails. At $U=7$, $10Q\equiv k\bmod7$ would put the nonsquare $3Q$ among the squares 1 and 2.

For arbitrary exponents at ten, (1) leaves just seven pairs of residue classes. In the following table the entries in each triple are $(x,P\bmod 5U,Q\bmod 5U)$.

| $U$ | Residue triples |
|---:|---|
| 3 | $(3,14,13),(3,7,11),(3,11,7),(6,14,4),(6,7,8)$ |
| 7 | $(7,34,12),(7,17,24)$ |

For $P,Q$ both coprime to 105 the rows are disjoint, since $xPQ\equiv1\bmod5$ distinguishes the three values of $x$. The five pairs modulo 15 contribute $5/8^2$ of the unit-class pairs and the two modulo 35 contribute $2/24^2$. Unitarity further requires $U^2\nmid 1+10PQ$, retaining exactly $(U-1)/U$ of the lifts in each case. Thus **at least $635/672$ of the ordered pairs of unit residue classes modulo $2205=5\cdot9\cdot49$ certify every nested $P,Q$ in those classes**. This proportion counts residue classes. The verifier independently enumerates all 1,016,064 pairs of unit classes.

Run **python3 check_nested_power_residues.py** to generate the table. A second calculation tries all scale divisors and actual powers in the combined modulus for every exclusion in the table and every $b\le40$. The finite congruence exclusions for the two cofactor-ten families also pass Lean 4.36. [Python verifier](check_nested_power_residues.py) · [Lean residue certificates](NestedPowerResidues.lean)
