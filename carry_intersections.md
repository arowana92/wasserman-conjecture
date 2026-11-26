# Simultaneous carry tests in different bases

Let $Q$ denote a prime power greater than 1 and let $1\le d\le1000$ with $\gcd(Q,d)=1$.

1. **The full conjecture holds at every $N=dQ+1$ with $Q\le10^{1000}$.** Over unbounded $Q$, the number of possible exceptions $N\le X$ in the union of all these families is $O_\epsilon(X^\epsilon)$ for each $\epsilon>0$.
2. **The full conjecture holds at every $N=5p^a q^b+1$ through $5\cdot10^{1000}+1$, and at every $N=7p^a q^b+1$ through $7\cdot10^{1000}+1$.** Here $p,q$ are primes and $a,b\ge1$. Over unbounded weights the numbers of possible exceptions up to $X$ in these two families are $O(X^{0.101})$ and $O(X^{0.212})$, respectively.

The [three-base argument](five_seven_pairs.md) strengthens both exceptional-set bounds in item 2 to $O_\epsilon(H^\epsilon)$ in every interval of length $H$.

The two finite bounds in item 2 concern the product $p^a q^b$ with both factors free to vary. The [elementary digit arguments](predecessor_five_seven.md) also give independent finite checks and an $O((\log X)^6)$ exceptional-set bound in the even factor-seven family. All conclusions use the [all-weight reduction](all_weights_endpoint_reduction.md).

## Necessary additions

For item 1, the earlier [cofactor sieve](prime_power_cofactors.md) excludes every even $Q$ in this range of $d$. For odd $Q$ it leaves 16 possible row shapes at 15 cofactors. A shape has
$$A=(aQ,bQ,cQ+1),\qquad C=V(f,g,h),\qquad UV=N=dQ+1,\qquad U=f+g+h.$$
Every prime dividing $U$ must be avoided by $A$. At a full power $w$ in $d$, a row avoiding its base has two parts divisible by $w$. Since $\gcd(V,w)=1$, if fewer than two of $f,g,h$ are multiples of $w$, its base must also be avoided by $A$. The checker regenerates the surviving shapes and uses all these forced bases. It imposes $dQ\equiv-1\bmod U$.

For item 2 write $t=(N-1)/5$ or $(N-1)/7$. In the only remaining factor-five case $A$ has part $3t$ and its complement sums to $2t+1$. Adding these two integers must be carry-free in bases 2, 3 and 5, and $t\equiv7\bmod15$. In the even factor-seven case the analogous summands are $6t,t+1$, the bases are 2, 3 and 7, and $t\equiv11\bmod72$. These tests follow from the [factor-five and factor-seven proofs](predecessor_five_seven.md).

At an odd factor-seven weight, $U$ is odd, the full factor 7 gives $U\mid k<7$, and the associated part $kt$ is 1 modulo 7. Hence $(U,k)$ equals $(3,3),(3,6)$ or $(5,5)$. The two summands $kt,(7-k)t+1$ must add without carries in bases $2,U,7$, since $A$ takes all bases of $N-1$ and the prime $U$. Also $t$ is even and $7t+1\equiv0\bmod U$.

## Exhaustive search by intervals

For fixed nonnegative slopes and offsets, testing whether entries $a_i t+b_i$ add without carries in base $p$ is a finite-state computation. Initialize the state at the offsets. On input digit $v$ read from least to most significant, compute $s_i+a_i v$ at each position. If the residues modulo $p$ sum to at least $p$, discard this transition. Otherwise the successor is $(\lfloor(s_i+a_i v)/p\rfloor)_i$. Once the input digits are exhausted, accept exactly if the components of the remaining state add without carries. This follows by induction from multiplication with carry. Each component of a state stays bounded by $\max(a_i,b_i)$.

To intersect several bases exactly up to a bound $B$, pad each input in each base to a fixed length with $p^m>B$. A high-digit prefix in that base specifies an integer interval $[j p^r,(j+1)p^r-1]$. Read such a prefix backward through the transitions, beginning with the set of accepting states. Keep only states reachable from the initial state in exactly $r$ unrestricted low digits. An empty set proves that no input in that interval can pass this base.

Start with $[2,B]$ and empty prefixes. In one base at a time, subdivide the widest current prefix interval by its next digit. Intersect its children with the interval already retained by all the bases. Each integer in the parent belongs to exactly one child in this subdivision. Discard empty intersections and impossible states. When at most 32 integers remain, test the congruence and all the actual additions directly. Thus every valid input is either returned or covered by an explicit failed condition.

For $B=10^{1000}$ the 16 cofactor shapes used 596,654 interval nodes and 89 direct checks in total. The factor-five test used 118,280 nodes and 31 direct checks. The even factor-seven test used 496,488 nodes and 26 direct checks. All three searches returned an empty set. The small-input cross-check compares with exact multinomial products and includes accepted inputs.

At odd factor-seven weights only $t=8452$ passes its two-summand tests, with $(U,k)=(5,5)$. Here $N-1=4\cdot7\cdot2113$. The unique three CRT parts are $44373,42260,31696$. Their sum is $118329$, whereas $N=59165$. Thus they cannot give the required row.

## Counting at unbounded weights

For a carry graph in base $p$, let $T_p(m)$ be the count of paths of length $m$ from all reachable starting states, with any ending state. Let the allowed digits in base $p^m$ be the values of all readable blocks of $m$ digits. Every admissible input belongs to the resulting integer missing-digit set. The number of different allowed digits is at most $T_p(m)$. To justify the containment, partition the padded expansion of the input into blocks. Each block occurs along an actual path, and the integer value of that block has its usual place-value order. Padding at the top with zero digits is valid for every accepted input. The all-zero block is allowed in these graphs. The integer missing-digit set has dimension at most $\log T_p(m)/(m\log p)$, by counting strings of blocks.

Theorem B of Glasscock, Moreira and Richter applies to integer missing-digit sets in multiplicatively independent bases. If their dimensions are at most $\alpha,\beta$, the intersection has at most $O_\eta(X^{\max(0,\alpha+\beta-1)+\eta})$ elements through $X$ for each $\eta>0$. The bases $p^m,q^m$ for distinct primes $p,q$ satisfy the independence hypothesis. [Glasscock, Moreira and Richter, *Additive and geometric transversality of fractal sets in the integers*](https://arxiv.org/abs/2007.05480)

At each of the 16 cofactor shapes the checker finds two forced bases with dimension bounds summing to at most $86/100$. It uses $m=64$ and bounds a dimension by $j/100$ only after checking the integer inequality $T_p(m)^{100}\le p^{mj}$. The two-base theorem gives $O_\epsilon(X^\epsilon)$ per shape, hence for their finite union.

For the factor-five family use bases 2 and 5 on the same parameter $t$, with slopes $(3,2)$ and offsets $(0,1)$. At $m=16384$ the path counts have bit lengths 9789 and 19139 respectively. Since $5^{1000}>2^{2321}$, the sum of the dimension bounds is strictly below
$$\frac{9789+(1000/2321)19139}{16384}<1.101.$$
The strict margin absorbs the $\eta$ in Theorem B, giving the asserted $O(X^{0.101})$ bound.

For the odd factor-seven cases with $k=3,5$, use bases 2 and 7 in the same calculation. The pairs of bit lengths are $(9039,28087)$ and $(10331,26697)$ respectively. Using $7^{1000}>2^{2807}$ gives dimension sums strictly below $1.212$ in both cases. For $k=6$ the monotone base-seven proof applies to $t$ at both parities, since $6t\equiv1\bmod7$ and $A$ avoids 7. Its count is $O((\log X)^6)$. This proves the second exceptional-set bound.

Run `python3 check_carry_intersections.py` from the project directory. The parameter `--digits L` sets the upper limit to $10^L$. The default is $L=1000$. The earlier sieve covering every cofactor through 1000 is reproducible with `python3 check_prime_power_cofactors.py 1000`.
