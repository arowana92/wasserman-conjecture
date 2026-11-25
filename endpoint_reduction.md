# Exact search when one row is (1,1,N−2)

This section concerns the selection containing the specified row. Its coefficient is $N(N-1)$. The other two coefficients form a counterexample with this row if, between them, they avoid every prime dividing $N(N-1)$. In particular, no factorization of a full trinomial coefficient is needed to enumerate candidates.

## Splitting the two consecutive integers

For each prime dividing $N(N-1)$ assign its full power in the corresponding integer to one of the two rows that avoids the prime. This partitions the prime powers into coprime products $UV=N$ and $ST=N-1$, with $U,S$ assigned to the first row and $V,T$ to the second. If both rows avoid a prime, either assignment is allowed.

A first row avoiding the primes assigned from $N$ has all parts divisible by $U$, so it has the form $(Ux_1,Ux_2,Ux_3)$ with $\sum x_i=V$. The second row has parts $Vy_i$ with $\sum y_i=U$. In particular $U,V\ge3$.

For each prime power $w$ assigned from $N-1$ to the first row, avoidance requires the residues of the actual parts modulo $w$ to be $1,0,0$ in some order. Since $UV\equiv1\bmod w$, the normalized parts have residues $V,0,0$ modulo $w$. Assign the exceptional residue to a bin for each $w\mid S$ in the factorization. The Chinese remainder theorem determines $x_1,x_2,x_3$ modulo $S$. Every possible row is obtained by taking the least positive representatives and adding nonnegative multiples of $S$ that make their sum $V$. The case $S=1$ gives all positive compositions of $V$. Interchange $U,V$ and $S,T$ for the other row. The program checks the full base-$p$ addition after generating these candidates, since the low-digit conditions alone are insufficient.

At least one of $S\ge V$ or $T\ge U$ holds. Otherwise $ST\le(V-1)(U-1)<N-1$. Both cannot hold since $UV=N>ST$. Name the first row so that $S\ge V$. Each bin must receive an exceptional residue for some prime power of $S$, since a bin receiving none would make its positive $x_i$ a multiple of $S$ even though $x_i<V$. Thus the first row must be assigned **at least three distinct** prime factors of $N-1$, with all three exceptional positions occurring. All its normalized parts are determined by their residues, without any lifting. In particular $\omega(N-1)\le2$ excludes the specified row, in agreement with Bergman’s Proposition 7. [Bergman](https://arxiv.org/abs/0806.0607)

Two further necessary inequalities prune the search. For each assigned prime power $w$ in $S$, two of the $x_i$ are multiples of $w$, so $2w+1\le V$. Also $S^2\mid x_1x_2x_3\le(V/3)^3$. Both arguments apply to the second row with the corresponding factors.

**Finite result.** The exhaustive search found no counterexample containing $(1,1,N-2)$ for any $3\le N\le1{,}000{,}000$. The independent small-input check matched all 1,564 assignments through 120, including 425 feasible assignments.

## Reproduction

Run `python3 search_endpoint.py 3 1000000` with Python 3.9 or later. For a surviving pair the program prints all three rows and checks their multinomial gcd using exact binomial products. The search covers only selections containing $(1,1,N-2)$. The [weighted obstruction](weighted_obstruction.md) gives a sufficient prime-gap hypothesis for this search to cover an entire weight.

The independent small-input check `python3 check_endpoint_search.py` tests each factor assignment separately, including feasible assignments. It compares the CRT search with exact coefficient enumeration before combining the two rows.
