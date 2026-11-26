# Sparse odd bases and binary carries

Two primes are *incompatible* at $N$ if every positive trinomial coefficient of weight $N$ is divisible by at least one of them. Write $s_q(N)$ for the base-$q$ digit sum.

**At an even weight, 2 is incompatible with every odd prime $q$ such that $s_q(N)\le6$.** Divisibility of $N$ by $q$ is optional. For $q\equiv1\bmod4$ the threshold rises to eight. These bounds are sharp even when $q\mid N$. The row $(12,162,768)$ at 942 avoids 2 and 3 with $s_3(942)=8$, and $(80,130,260)$ at 470 avoids 2 and 5 with $s_5(470)=10$.

**If $4\mid N$ and two distinct odd primes $p,q\mid N$ have $s_p(N),s_q(N)\le6$, the three primes $2,p,q$ certify every pair of rows together with the endpoint.** The [endpoint reduction](all_weights_endpoint_reduction.md) proves the full conjecture at every such weight. The same three-prime certificate works at any even weight under either of the following conditions on these odd primes. One of the digit sums is at most five, or at least one of the orders of $p\bmod q$ and $q\bmod p$ is odd. The latter includes any two primes $p,q\equiv3\bmod4$.

**If $N$ is even, $s_2(N)\le6$, and $s_q(N)\le6$ for some odd prime $q\mid N$, the full conjecture holds.** Indeed the [single-coefficient theorem](six_units_weight_prime.md) forces any counterexample to have one row avoiding both sparse weight primes, contradicting the binary incompatibility. Together with the [mixed theorem](mixed_six_units.md), at an even weight of binary digit sum at most six it suffices that any odd prime of $N(N-1)$ has digit sum at most six.

## Counting the units of even parts

In a binary carry-free row, the least occupied bit of each part has a different position. The number of odd parts is $\eta=N\bmod2$. If the row also avoids an odd prime $q$, its $s_q(N)$ units are split among the parts. A positive even part requires an even number of units. With just two it has the form $q^a(1+q^d)$ for $a,d\ge0$, and its 2-adic valuation is either 1 or $h=v_2(q+1)$. This follows by reduction modulo 4 for even $d$ and odd-power factorization for odd $d$.

Let $J(N)$ be the set of occupied binary positions of $N$, numbering the units place zero. At most $|J(N)\cap\{1,h\}|$ even parts can therefore use just two units. Each other even part uses at least four and each odd part at least one. Thus simultaneous avoidance requires
$$
s_q(N)\ \ge\ \eta+4(3-\eta)-2\,|J(N)\cap\{1,h\}|.
$$
For even $N$ this gives the claimed bounds. If $4\mid N$, it also gives $s_q(N)\ge10$ at every odd base, and at least 12 if the bit at $h$ is absent.

A complementary bound uses $c=v_2(q-1)$. If $e_1<e_2<e_3$ are the three lowest occupied binary positions of $N$, then simultaneous avoidance requires
$$s_q(N)\ \ge\ \sum_{i=1}^3 2^{\min(e_i,c)}.$$
Indeed $s_q(z)\equiv z\bmod2^c$, so a positive part $z$ costs at least $2^{\min(v_2(z),c)}$ units. Its ordered valuations are at least the respective $e_i$. With fewer than three occupied positions all positive rows already carry at 2.

## Two odd prime factors of the weight

A row avoiding $p,q\mid N$ has every part divisible by both. Thus no part is a singleton power in either base, and at most six units require exactly two for each part in each base. Each part has the form
$$
p^{a_i}(1+p^{d_i})=q^{b_i}(1+q^{f_i}),\qquad a_i,b_i,d_i,f_i>0 .
$$
The gaps are positive since $p,q$ are distinct and odd. As $p^{d_i}\equiv-1\bmod q$ and $q^{f_i}\equiv-1\bmod p$, both orders must be even. Each $d_i$ is an odd multiple of half the order of $p$ modulo $q$. All three parts consequently have the same 2-adic valuation, equal to $v_2(N)$ since their sum has three odd units after scaling. If $4\mid N$, this requires $p\equiv q\equiv3\bmod4$. At two such primes the orders cannot both be even. An even order at a prime congruent to 3 modulo 4 makes the quadratic character $-1$, and quadratic reciprocity requires opposite characters. This proves the pairwise incompatibilities. Three pairwise incompatible primes dividing the endpoint coefficient certify the endpoint and every pair of other rows, since two rows can avoid at most two of the primes between them.

## The common cofactor at the remaining even weights

Suppose two weight primes $p,q$ have at most six units and a counterexample still exists. The case of an odd weight follows from the five-unit result, since odd bases have odd digit sums there. The restrictions above leave $N\equiv2\bmod4$ and odd $p,q$ jointly avoided by the same row $A$. Let $P,Q$ be their full powers in $N$. Put $g=\gcd_i d_i$, $h=\gcd_i f_i$. The gcd of the three parts satisfies
$$
G=P(1+p^g)=Q(1+q^h)=DPQ,\qquad D\equiv2\bmod4,\qquad \gcd(D,pq)=1 .
$$
Indeed $\min a_i=v_p(N)$ by carry freedom. The $d_i/g$ are all odd and have gcd one. Modulo the gcd of the three $1+p^{d_i}$, Bézout on these quotient exponents gives $p^g\equiv-1$ since the sum of the Bézout coefficients is odd. Conversely $1+p^g$ divides each factor. The same argument at $q$ gives the formula. The assigned full powers divide $G$, whose valuations at $p,q,2$ agree with those of $N$.

In fact **$D\ge6$**. Write $M=N-1$, with $U,S$ the products of weight and predecessor powers assigned to $A$ and $T=M/S$. The single-coefficient theorem gives $2US<M$ because $US$ divides $NM/\gcd(NM,K_A)$. The product of the three normalized companion parts gives $27T^2\le U^3$. Thus $U>108$, and $U\mid G$. If $D=2$, the equalities $1+p^g=2Q$, $1+q^h=2P$ are exactly the [reciprocal power equations](mixed_six_units.md#a-reciprocal-prime-power-calculation). Their only solution in either order gives $G=30$, a contradiction.

In particular **$N/(2PQ)$ must be composite**. It equals $(D/2)(N/G)$, with $D/2\ge3$ and $N/G$ the sum of three positive odd integers. Thus the two sparse prime factors also suffice whenever $N/(2PQ)$ is one or a prime.

## A minimum certificate at 396

At $396=4\cdot9\cdot11$, the primes $2,3,11$ are pairwise incompatible by the theorem. Both odd-base digit sums are six. A row avoiding 389 has two parts of sum at most 7. For binary avoidance the two smallest parts would sum to at least $4+8=12$. At 3 or 11 all parts would be divisible by 9 or 11, respectively. Hence $2,3,11,389$ are pairwise incompatible, giving a four-prime certificate directly for any three rows.

It is minimum. The three weight primes have digit sums greater than two, and a prime outside the weight having digit sum at most two would require $395=5\cdot79$ to be a power of that prime. Thus every individual prime is avoidable, by splitting at least three digit units among the parts. Three rows can defeat any list of at most three such primes.

The verifier checks 31,439 joint avoidances among 9,121 exact odd coefficients, all 13,068 unordered positive rows at 396, and 192,024 digit allocations for incompatible weight-prime pairs through 200,000. Run `python3 check_binary_sparse_incompatibility.py`.
