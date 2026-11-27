# Weights whose predecessor divides an initial least common multiple

Put $L_B=\operatorname{lcm}(1,\ldots,B)$. A positive integer divides $L_B$ exactly when every full prime power in its factorization is at most $B$. Write $M=N-1$ and let $K_A$ be the multinomial coefficient of a positive row $A$ of weight $N$.

- **If $M\mid L_{61}$ and $N\ge3$, Wasserman's three-row conjecture holds at $N$.** This covers all 3,538,943 such weights, the largest being $591133442051411133755680801$. More precisely, whenever $N=UV$, $U,V\ge3$, $\gcd(U,V)=1$, every pair of positive trinomial rows scaled by $U,V$ respectively has a common coefficient prime dividing $M$.
- **If $M\mid L_{53}$ and $N\ge3$, every individual coefficient from at least three positive parts satisfies $\gcd(M,K_A)\ge\gcd_i A_i$.** In particular $\gcd(NM,K_A)\ge N$. This covers all 884,735 such weights, up to $164249358725037825439201$.
- **For every $d\mid L_{61}$ and every prime $q$, all weights $1+dq^e$ satisfy the conjecture whenever $e$ is a positive multiple of the exponent of the unit group modulo $L_{d+1}$.** The simpler condition $(d+1)!\mid e$ suffices. If $d\mid L_{53}$, every $q$-avoiding coefficient from at least three positive parts at these weights satisfies $\gcd(d,K_A)\ge\gcd_i A_i$.

The finite assertions use exhaustive computation. The last assertion allows arbitrarily large base primes and exponents. For example one can take the same period $(L_{61}+1)!$ for all of the indicated cofactors and primes.

## Reducing the pair to one row

Suppose the first assertion about a pair fails. Write the two rows as $U x$ and $V y$, with normalized sums $V,U$. Assign each full power $p^e\parallel M$ to a row whose coefficient avoids $p$. Let $S,T$ be the two assigned products. One row must take a product at least as large as its normalized sum, since otherwise $M=ST\le(V-1)(U-1)<N-1$. Interchange names so that $S\ge V$. It therefore suffices to exclude a single such row at each unitary factorization. No carry tests at the primes of $N$ enter the computation.

To obtain the full conjecture, use the [endpoint reduction](all_weights_endpoint_reduction.md). In a failure, assign the full powers of $N$ to the two nonendpoint rows avoiding the respective primes. A row avoiding $p\mid N$ has every part divisible by the full power of $p$ in $N$. The assigned products $U,V$ are coprime, $UV=N$, and each is at least three since it is the normalized sum of the other row. The preceding result supplies the required common prime of $M$.

For the individual assertion one can take $U=\gcd_i A_i$ and $V=N/U$. The case $U=1$ is automatic. Suppose $J=\gcd(M,K_A)<U$. Then $M/J\ge V$, since $M=UV-1$. At each $p^e\parallel M$, retain the power $p^b$ with $b=e-\min(e,v_p(K_A))$. When $b>0$, all but one of the actual entries are zero modulo $p^b$ and the other is one. Indeed, an overflow in the sum of the residues at layer $j\le e$ persists through layer $e$ because the weight residue stays one. Each such layer adds at least one to the factorial valuation. Thus the product of the retained powers is at least $V$ and $v_p(K_A)\le e-b$ at every retained power. This time the computation uses **all** divisors $U\ge2$ of $N$ with $V\ge3$.

For any row of length greater than three, merge its parts into three nonempty groups. Its coefficient is a multiple of the grouped coefficient, and its part gcd divides the grouped part gcd. This proves the individual inequality for every indicated length. The decrement identity gives $N/U\mid K_A$, so the two coprime factors $\gcd(N,K_A)$ and $J$ have product at least $N$. For two rows with this bound, their individual endpoint gcds multiply to at least $N^2$, while their lcm divides $N(N-1)$. Hence they share a prime dividing the endpoint. [Lean verifies the divisor step](PrimitiveRowScales.lean).

### An individual bound with no prime-power ceiling

Suppose $A$ has $r$ positive parts, $U=\gcd_i A_i$, $V=N/U$, and $D>1$ divides $N-1$. If $D$ has $h<r$ distinct prime divisors, then
$$
           \gcd(D,K_A)\ \ge\ \frac{(r-h)D}{V-h}.
$$
Retain from $D$ the powers as in the individual proof above and put $T=D/\gcd(D,K_A)$. Each retained power has just one exceptional position. At least $r-h$ normalized entries are divisible by all of $T$, so $V\ge(r-h)T+h$. This proves the bound even when some powers of $D$ are proper subpowers of those in $N-1$.

In particular, if $N-1$ itself has $h<r$ distinct prime factors, then $\gcd(N-1,K_A)\ge(r-h)U$, with strict inequality when $hU>1$. This follows from $hU\ge1$ and $N=UV$. Thus the individual trinomial result also holds with **two completely unrestricted prime powers** in the predecessor. More generally every such coefficient from $r>h$ parts satisfies $\gcd(N(N-1),K_A)\ge(r-h)N$. The transfer below consequently gives $\gcd(d,K_A)\ge\gcd_i A_i$ at every eligible exponent when $d$ has at most two prime divisors, even outside $L_{53}$.

## Searching the assigned powers

Order the full powers of $M$ in decreasing order. For a hypothetical row, take the prefix of its assigned or retained powers up to the first product $P\ge V$. Each power $w$ used divides two positive normalized parts, so $2w+1\le V$. Every position must be exceptional at some selected power, since otherwise its normalized part would be a positive multiple of $P$ strictly below $V\le P$. In particular at least three distinct primes are selected. The [native verifier](check_smooth_predecessors.cpp) tries every feasible prefix. In the individual test it tries every retained exponent from 1 to $e$ at a selected prime.

For each position assign a subset of the selected powers. Its normalized part is congruent to $V$ at precisely its assigned powers and zero at the others, because $UV=N\equiv1\pmod P$. These congruences recover the part uniquely in $[0,P)$. First enumerate subsets for one smallest part $x$ with $0<x\le V/3$. Its actual part $a=Ux$ must pass $v_p\binom Na\le e-b$ at every selected prime. In this binomial cut the lowest $b$ digits have no borrow, so they may be discarded. In the full-power test the allowance is zero. Split the remaining powers between the remaining positions, reconstruct $y,z$, require $x+y+z=V$, and check the valuation of the second factor in $K_A=\binom Na\binom{N-a}{Uy}$. Every genuine row would pass the two tests on at least one enumerated prefix.

## Consequences at unbounded weights

For a common scale of two arbitrary primitive trinomial directions, suppose every pairwise coordinate gcd within each direction divides $L_{61}$. Then the two coefficients and $N(N-1)$ share a prime at **every** scale. Indeed the three pairwise gcds of a primitive triple are pairwise coprime, and their product is $\gcd_{i<j}(x_i x_j)$. An exception would force $N-1$ to divide the lcm of the corresponding two products by the [primitive row criterion](primitive_row_scales.md). Then $N-1\mid L_{61}$, already covered. This applies in particular if the two rows each satisfy a homogeneous linear relation whose three coefficients are nonzero divisors, up to sign, of $L_{61}$. In a primitive triple the gcd of two entries divides the coefficient on the remaining entry, by reducing the relation modulo that gcd.

For the weights $N=1+dQ$, $Q=q^e$, use the [specialization to weight $n=1+d$](small_weight_transfer.md#specializing-the-two-rows). The exponent hypothesis gives $Q\equiv1\pmod{p^j}$ for every $p^j\le n$ with $p\ne q$. In a counterexample that specialization produces a unitary factorization $n=UV$ with $U,V\ge3$ and two rows scaled by $U,V$. Every prime of $d$ is avoided by at least one of the two small rows. Our finite theorem rules this out even if carry information at all weight primes is discarded. At $d=1$, the row avoiding $q$ would require two positive multiples of $Q$ among its parts, already impossible.

For the stronger individual conclusion consider three parts first. In a $q$-avoiding row they have the form $aQ,bQ,cQ+1$ with $a,b\ge1$ and $a+b+c=d$. Let $U$ be their gcd. If $U>1$ then $c\ge1$, $U\mid a,b$ since $q\nmid U$, and $U\le(d-1)/2$. Every full power of $U$ is at most $d+1$ and coprime to $q$, so $Q\equiv1\pmod U$. The small row $A_0=(a,b,c+1)$ has gcd a multiple of $U$. At each $p\mid d$ distinct from $q$, all valuation contributions from layers $p^j\le n$ agree in the two rows by the congruences on $Q$. Larger layers contribute zero at $n$ and nonnegatively at $N$. If $q\mid d$, the digits of $a,b,c$ add carry-free to $d$, since they are precisely the digits above the lowest $e$ places of the large row. Their units digits are zero, so $a,b,c+1$ still add carry-free. Thus $\gcd(d,K_{A_0})\mid\gcd(d,K_A)$. The individual inequality at $n$ proves the claim. Grouping proves it at every length of at least three. In particular $\gcd(Nd,K_A)\ge N$ for those $q$-avoiding coefficients.

## Reproduction and scope

The [input generator](check_smooth_predecessors.py) enumerates every divisor of $L_B$ except 1 and checks the primality and product of the factors used. The [independent verifier](check_smooth_predecessors_cross.py) reconstructs actual parts from modular inverses of their divisors and uses Legendre factorial floors. On 2,011 weights it reproduced all 27,884 full-power prefixes, 25 first cuts and two rows, and all 68,194 retained-power prefixes, 56 first cuts and two rows. It also checked 114,831 small exact multinomials, 9,532 bounds using divisors of the predecessor, and 1,029 transferred rows, including 90 with nontrivial part gcd.

The positive control at $N=222788354$ records a genuine limitation on the individual inequality. For $A=(45491401,86617810,90679143)$ the part gcd is 283, while $\gcd(N-1,K_A)=251$ and $\gcd(N(N-1),K_A)=197596738<N$. Here $N-1=29\cdot127\cdot241\cdot251$, outside both divisor families. The full conjecture at this weight is covered by the [initial range computation](full_range_1b.md).

With python-flint installed, run from the repository root using a C++17 compiler with unsigned 128-bit integer support.

~~~~sh
mkdir -p .work
c++ -O3 -std=c++17 check_smooth_predecessors.cpp -o .work/smooth
python3 check_smooth_predecessors_cross.py .work/smooth
python3 check_smooth_predecessors.py .work/smooth --bound 53 --minimum-predecessor 2 --retained --jobs 6
python3 check_smooth_predecessors.py .work/smooth --bound 61 --minimum-predecessor 2 --jobs 8
~~~~

| Computation | Weights | Prefixes | Smallest cuts | Whole rows |
|---|---:|---:|---:|---:|
| Full powers through 61 | 3,538,943 | 3,488,544,028 | 1,580 | 0 |
| Retained powers through 53 | 884,735 | 1,520,290,124 | 1,679 | 0 |

The input SHA-256 digests in the indicated order are `5fc1346cfcfdb99579c869f0882ae08f6a530bc49ea5e31cc88288d88a8984c6` and `ec09cd714204ac82d713268379052889f39a540c30fcfd198116b3a655205dc8`. The full-power calculation was also run separately above and below predecessor $10^9$. The two blocks returned 322 and 1,258 smallest cuts respectively, with no whole row in either block.
