# Which single cuts can fail against a trinomial?

For each $N\ge3$ let $W_N$ be the set of $1\le a\le N/2$ such that at least one positive trinomial row $B$ of weight $N$ satisfies
$$\gcd\left(N(N-1),\binom Na,K_B\right)=1.$$
Reflection identifies the cuts $a$ and $N-a$. Membership allows a different $B$ for each cut.

**Finite theorem.** Exactly 20 weights through $1{,}000{,}000$ have $W_N\ne\varnothing$. At each such weight there is only one cut up to reflection, as listed below. In particular, throughout this range **every positive trinomial row has a cut whose binomial coefficient shares an endpoint prime with every positive trinomial row of that weight**. The prime may depend on the second row. This is a strengthening of the usual two-row endpoint assertion in this range.

| $N$ | $a$ | $N$ | $a$ |
|---:|---:|---:|---:|
| 5,752 | 568 | 162,428 | 8,260 |
| 12,273 | 768 | 191,417 | 19,384 |
| 20,465 | 1,280 | 302,057 | 6,664 |
| 21,809 | 5,424 | 327,113 | 96,287 |
| 59,165 | 16,905 | 392,684 | 92,396 |
| 60,524 | 2,088 | 514,484 | 6,416 |
| 72,989 | 5,397 | 848,897 | 17,409 |
| 95,657 | 25,001 | 896,888 | 272,472 |
| 105,596 | 2,156 | 915,932 | 100,376 |
| 121,076 | 4,176 | 920,768 | 264,256 |

Indeed if all three parts of a row lie in $\{a,N-a\}$ and sum to $N$, they must all equal $N/3$. None of the listed cuts has $3a=N$. Any part outside that pair furnishes the asserted cut.

## The first failure and its refinements

At $N=5752=8\cdot719$, the endpoint coefficient is $N(N-1)=8\cdot719\cdot81\cdot71$. For $a=568$ and $B=(719,719,4314)$, the respective gcds of the binomial and trinomial coefficients with the endpoint coefficient are **719** and **46,008**. They are coprime. The independent exact search finds precisely five witnessing rows at this first cut, the multiples of 719 of the five unordered positive triples of sum eight.

Every refinement of either side of the cut at 568 restores at least two common endpoint primes with the displayed $B$. The counts for common-prime sets $\{2,3\}$, $\{2,71\}$, $\{3,71\}$, and $\{2,3,71\}$ are respectively 9, 10, 10, and 2,847, counting each unordered split of the chosen side once.

The two-prime lower bound also follows from the digits. For the side 568, cuts in the lower half avoiding 2 are nonzero subset sums of $8,16,32$. Those avoiding 3 belong to $\{1,81,82,243,244\}$, since $568=2\cdot243+81+1$. Those avoiding 71 are multiples of 71. The three sets are pairwise disjoint. For the side $5184=64\cdot81=71^2+2\cdot71+1$, cuts avoiding 2 or 3 must be multiples of 64 or 81 respectively. No interior cut is a multiple of both; the positive lower-half choices avoiding 71 are $1,71,72,142,143$, none divisible by 64 or 81. Kummer's criterion proves the claim for both sides.

## An all-weight bound with two predecessor primes

Fix $B$. Assign every endpoint prime avoided by $B$ to it and the rest to the cut. Write $UV=N$, $ST=M=N-1$, with $V,T$ the products of the full powers assigned to $B$ from the weight and predecessor. A witnessing cut has $a=Ux$ and $0<x<V$. Each full power $w$ in $S$ requires $x\equiv0$ or $V\pmod w$. At each power assigned to $B$, two normalized entries of $B/V$ are zero modulo that power.

If $T$ has at most two prime bases, one normalized entry of $B$ is divisible by all of $T$, since at most two positions receive the exceptional residues. Thus $U>T$, whence $S=(UV-1)/T>V$ for a witnessing cut ($V\ge2$). Each choice of a subset of the primes of $S$ now determines at most one $x$ in $(0,V)$. The empty and full subsets give none. Reflection exchanges complementary subsets. For $s=\omega(S)\ge1$ this gives at most $2^{s-1}-1$ unoriented cuts for this $B$; when $S=1$ the premise $S>V$ is impossible.

**For every weight with at most two distinct primes dividing $M$, $|W_N|\le1$ even as $B$ varies.** A witness needs at least two primes in $S$ by the preceding result. Thus $S=M$ and both predecessor primes are assigned to the cut. The original index $a$ must select a nonempty proper subset of these two powers at which it is one, being zero at the other. Modulo $M$ those two CRT residues sum to $1\bmod M$ and give a single complementary pair in $(1,N-1)$. The boundary cut 1 cannot witness: a trinomial coefficient coprime to $N$ would have all its positive parts divisible by $N$, using $N\mid b_i K_B$. Finally the exceptional cut cannot equal $N/3$, since $M\mid a(a-1)$ and $3a=N$ would imply $M\mid2$. This proves the stronger cut assertion at all these weights.

## A criterion equivalent to the original problem

For a fixed $B$ take all primes dividing both $K_B$ and $NM$. Say $x\preceq_B y$ if every digit of $x$ is at most the corresponding digit of $y$ in each of those bases. The positive interior elements below $N$ are precisely the oriented cuts that witness against $B$. A counterexample containing the endpoint and $B$ exists **if and only if two distinct such cuts are comparable**. Indeed $0<x<y<N$ with $x\preceq_B y\preceq_B N$ gives the carry-free row $(x,y-x,N-y)$ simultaneously at every required prime. Conversely the partial sums of a row avoiding all the primes give just such a chain. By the [endpoint reduction](all_weights_endpoint_reduction.md) this covers the general conjecture. The unrestricted comparability question remains open here.

## Reproduction

The main verifier splits the full powers in both $N$ and $M$ between the cut and the trinomial. For the cut it enumerates all its CRT choices modulo $S$ and all lifts with $1\le x\le V/2$, then checks all assigned bases. Since $UV\equiv1\bmod S$, the two possible residues of $x$ at each power are zero and $V$. For each assignment with any retained cuts, the [trinomial verifier](search_endpoint.py) decides whether a row exists. A witness for the assignment certifies all its retained cuts. The early bounds $2w<U$ and $27T^2\le U^3$ follow respectively from two positive multiples of $w$ in $B/V$ and from $T^2$ dividing its product of parts. A single weight-prime factor forces a carry in every interior binomial and a single predecessor-prime factor allows only the boundary cuts in the binomial if the trinomial carries there, as it must. The program checks an explicit trinomial witness for every retained cut by exact coefficients, reports every $W_N$, and tests whether its oriented cuts can sum to $N$ in triples.

The second verifier independently enumerates all cuts and compatible positive rows through the first failure using exact binomial recurrences, products and gcds. It uses no prime factorization, digit test or residue assignment. At $5752$ it recovers all five rows.

~~~sh
python3 check_cut_relaxation.py        # uses python-flint
python3 check_cut_relaxation_exact.py  # standard library
~~~
