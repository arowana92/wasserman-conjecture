# The conjecture holds at almost every weight

**Theorem.** Let $E(X)$ count the integers $3\le N\le X$ for which some three positive trinomial coefficients of weight $N$ have gcd one. Then $E(X)=O(X^{999/1000})$. In particular the weights satisfying Wasserman’s full three-row conjecture have natural density one.

A [refinement](endpoint_count.md) improves the exponent to $15/16$.

## How a single part restricts the weight

First suppose the selection contains the endpoint $(1,1,N-2)$. Its coefficient $N(N-1)$ is even, so one of the other rows, say $B$, must be carry-free in base 2. Assign each prime of $N(N-1)$ to one of the two non-endpoint rows avoiding it. Let $U,S$ be the products of the full powers assigned to $B$ from $N,N-1$, respectively. Let $V,T$ be the complementary products, so $UV=N$ and $ST=N-1$. If both rows avoid a prime, choose either. [Endpoint assignment](endpoint_reduction.md)

The parts of $B$ are multiples of $U$ and those of the other row are $V y_1,V y_2,V y_3$ with $y_1+y_2+y_3=U$. Hence $U,V\ge3$. At a full prime power in $T$, two of the other row’s parts must be 0 modulo this power, since the residues add to 1. The same holds for the $y_i$ since $V$ is coprime to $T$. Thus $T^2\mid y_1y_2y_3\le(U/3)^3$, which implies
$$(US)^3\ge U^3 S^2\ge27(N-1)^2.\tag{1}$$

For **every** part $z$ of $B$ we have $U\mid z$ and $S\mid z(z-1)$. The latter follows because the three residues at each full power in $S$ are $1,0,0$. Here $z\ge U\ge3$, so both divisor conditions concern positive integers. Given $z\le 2X$, there are $X^{o(1)}$ possible pairs $(U,S)$, by the elementary divisor bound $\tau(n)\ll_\epsilon n^\epsilon$. This bound follows by using $e+1\le p^{\epsilon e}$ for all sufficiently large primes $p$ and absorbing a bounded factor for each remaining prime.

For each pair, the congruences $N\equiv0\pmod U$ and $N\equiv1\pmod S$ specify one residue class modulo $US$, since $U,S$ are coprime. For $X\le N<2X$, (1) bounds the number in that class by $O(X^{1/3})$. Thus each possible value of a single part $z$ accounts for at most $X^{1/3+o(1)}$ weights.

## Counting the binary digits

For convenience use blocks $2^{L-1}\le N<2^L$. Write $s_2$ for the number of binary ones and $H(t)=-t\log_2t-(1-t)\log_2(1-t)$ for binary entropy. The number of length-$L$ bit strings with at most $tL$ ones, for $0<t<1/2$, is at most $2^{L H(t)}$. This follows by bounding their total binomial probability at parameter $t$. By complementing bits the corresponding upper-tail bound holds for $t>1/2$.

The weights with $s_2(N)>13L/25$ number at most $2^{L H(13/25)}$. At every remaining endpoint counterexample, the three parts of $B$ use disjoint binary digits. Their digit counts sum to $s_2(N)$, so one part $z$ has $s_2(z)\le13L/75$. There are at most $2^{L H(13/75)}$ choices for it. The modulus count bounds the number of endpoint counterexample weights in the block by
$$2^{L H(13/25)}+2^{L(H(13/75)+1/3+o(1))}.$$
Both exponents are strictly below $999/1000$. The following exact integer inequalities verify the margins without rounding logarithms.
$$25^{1000}<2^{999}13^{520}12^{480},\qquad
75^{3000}<2^{1997}13^{520}62^{2480}.$$

## Selections without the endpoint

The [second-part bound](power_intervals.md) gives $b>N^{1/3}$ for the second-smallest part of every row of a sufficiently large counterexample without the endpoint. Let $p_n\le N<p_{n+1}$ be consecutive primes. A row avoiding $p_n$ has a part at least $p_n$ since $p_n>N/2$. Its two smaller parts therefore sum to at most $N-p_n$. Thus every such weight in $[X,2X)$ lies in a prime gap of length $d_n=p_{n+1}-p_n>X^{1/3}$.

The Peck–Maynard mean-square gap bound gives $\sum_{p_n\le2X}d_n^2\ll_\epsilon X^{5/4+\epsilon}$. The number of integers in the relevant gaps is at most $X^{-1/3}\sum_{p_n\le2X}d_n^2\ll_\epsilon X^{11/12+\epsilon}$. This includes a gap crossing either boundary of the block. [Maynard’s paper and attribution to Peck](https://arxiv.org/abs/1201.1787)

Summing the endpoint and non-endpoint bounds over dyadic blocks proves the theorem. The proof counts exceptional weights and leaves the existence of an individual exceptional weight open.
