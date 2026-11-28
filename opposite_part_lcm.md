# A restriction on the least common multiple of parts in different rows

Take two positive trinomial rows of weight $N$. Suppose that for every full prime power $p^e\parallel N(N-1)$ at least one of the rows has sum of its three residues modulo $p^e$ equal to $N\bmod p^e$. The assumption holds if the two coefficients and the endpoint coefficient $N(N-1)$ have gcd one. It uses only the lowest prime-power layers, so the theorem applies even when there is a common carry at a higher layer.

Choose any **single part** $a$ of the first row and $b$ of the second. Then both numbers
$$m=\frac{\operatorname{lcm}(a,b)}N,\qquad
k=m\frac{(N-a)(N-b)}{N(N-1)}$$
are positive integers. There exists a factorization $m=ABCD$ into four pairwise coprime positive integers (allowing factors 1), with
$$(A+B)(C+D)\le k<m. \tag{1}$$
Equivalently the minimum of the left side over all distributions of the full prime powers of $m$ into four groups is a computable lower bound on $k$.

In particular $16m\le k^2$ and $2(P+1)\le k$ for each full prime power $P\parallel m$. The quotient $m$ is at least 20 and is neither a prime power nor twice a prime power. If $m=uv$ consists of exactly two full prime powers, then $k\ge2(u+v)$. At $m=20$, for example, $k\ge18$ requires both selected parts to be smaller than $N/10$ when $a,b>1$. These conditions apply simultaneously to all nine pairs of single parts in a counterexample, without any upper bound on $N$.

## Proof

Assign each prime of $N(N-1)$ to a row passing its full power. Let $U,V$ be the products of the powers of $N$ assigned to the first and second row. Then $UV=N$, $U$ divides every part of the first row and $V$ every part of the second. In particular $N\mid\operatorname{lcm}(a,b)$.

At a predecessor power, the residues in its assigned row are a permutation of $(1,0,0)$. Let $S_1,S_2$ be the products of the powers assigned to the first or second row respectively at which its selected part has residue 1. Let $T$ contain the remaining powers, so $S_1S_2T=N-1$. Each power in $T$ divides $a$ or $b$ and is coprime to $N$, hence $T\mid m$. Write $d=m/T$.

The two companions of $a$ are positive multiples of $US_1$, and those of $b$ are positive multiples of $VS_2$. Let their normalized pairs be $i_1,i_2$ and $j_1,j_2$, with sums $u,v$. Then
$$uv=\frac{(N-a)(N-b)}{N S_1S_2}
=T\frac{(N-a)(N-b)}{N(N-1)},\qquad k=duv<m.$$
At each prime power of $T$, its assigned row's selected part has residue zero, so one of its two companions also has residue zero. All four normalizing factors $U,V,S_1,S_2$ are coprime to $T$. Assign the power to one such normalized companion. This yields four pairwise coprime products $t_1t_2t_3t_4=T$ dividing $i_1,i_2,j_1,j_2$, respectively. Consequently
$$\frac{k}{m}=\frac{uv}{T}\ \ge\
(1/t_1+1/t_2)(1/t_3+1/t_4).$$
Enlarge each assigned prime power to its full power in $m$, and distribute the full powers of any remaining primes of $m$ among the four products. This can only decrease the last expression. Multiplication by $m$ proves (1).

The arithmetic–geometric mean inequality gives $(A+B)(C+D)\ge4\sqrt m$. The group containing $P$ has size at least $P$, its partner at least 1, and the other pair sums to at least 2. For two full powers $u,v$, enumeration of the placements gives the minimum $2(u+v)$. These bounds exclude a prime power or twice a prime power. They exclude every $m\le19$, since $16m\le k^2<m^2$ requires $m>16$ and $17,18,19$ have already been excluded. For $a,b>1$ the fraction $k/m$ is strictly less than both $(N-a)/N$ and $(N-b)/N$.

The integers $m,k$ here coincide with the parameters in the [cross-cut theorem](cross_cut_ratios.md) when $D=J=1$. In particular this restriction on the two single parts applies at every product $hL$ considered in the [ratio enumeration](cross_cut_critical.md).

## Numerical check

Run `python3 check_opposite_part_lcm.py 10000`. For each weight through 10,000 with at least three distinct predecessor primes, the program allocates the full prime powers of $N$ and $N-1$, solves the three residue conditions in each row, and selects a positive representative of each feasible ordered residue pattern. It then checks the original full-power hypothesis on the actual rows and minimizes the four-factor cost independently over all multiplicative partitions for each of their nine part pairs. The resulting sample contains 41 row pairs at 35 weights and 369 choices of parts. Its least $m$ is 35, at $N=3256$ with rows $(280,1240,1736)$ and $(814,1221,1221)$.
