# One exponent progression for every prime and nearly a billion cofactors

Put $B_n=\operatorname{lcm}(1,\ldots,n)$ and let $L_n$ be the exponent of the group of units modulo $B_n$. **For every $1\le d<10^9$, every prime $q$, and every positive multiple $e$ of $L_{d+1}$, Wasserman's three-row conjecture holds at $1+dq^e$.** In particular, all positive multiples of the *same* integer $L_{10^9}$ work for every one of these cofactors and every prime. One may replace $L_n$ by $n!$ throughout.

There is no upper size restriction on $d$ in either of the following extensions.

- The theorem holds whenever $d$ has at most two distinct prime factors, or $d+1$ is a prime power or twice an odd prime power.
- Among $d\le X$, at most $O(X^{0.8689})$ can fail to satisfy the theorem for even one prime $q$ and one eligible exponent $e$.

The period can be enormous. The [earlier finite theorem](full_range_1b.md) concerns the weight itself. Here only the cofactor is bounded, and the congruence proof also gives every exponent for every prime $q\equiv1\pmod{B_{d+1}}$ when $d<10^9$.

## Specializing the two rows

Write $Q=q^e$, $N=1+dQ$ and $n=d+1$. It suffices that $Q\equiv1\pmod{p^k}$ for every prime power $p^k\le n$ with $p\ne q$. The exponent hypothesis gives precisely the congruences needed here, since every unit modulo $p^k$ lifts to a unit modulo $B_n$ by the Chinese remainder theorem. In particular $L_n\mid L_m$ for $m\ge n$, and $L_n\mid n!$ by Euler's theorem on each prime power.

In a counterexample the [endpoint theorem](all_weights_endpoint_reduction.md) gives the row $(1,1,N-2)$ with coefficient $N(N-1)$. Choose the other row $A$ that avoids $q$ and assign all full prime powers of $N$ to avoiding rows. The [factor assignment](shifted_directions.md#the-arithmetic-that-bounds-the-enumeration) gives
$$
 A=(aQ,bQ,cQ+1),\quad a+b+c=d,\quad U\mid a,b,\qquad
 C=V(f,g,h),\quad f+g+h=U,\quad UV=N.
$$
Both sets of multipliers are positive, $U,V\ge3$, $\gcd(U,V)=1$, and $U$ is assigned to $A$. Indeed the residues modulo $Q$ on $A$ are $0,0,1$ and a part equal to 1 would violate divisibility by $U$. Also $U\le(d-1)/2$ and $q\nmid U$. These facts already cover $d=1$.

Every prime power dividing $U$ divides $n$ by the congruences on $Q$. Moreover $U$ is a unitary divisor of $n$, meaning $\gcd(U,n/U)=1$. For $p^j$ exactly dividing $U$, we have $v_p(N)=j$. If $p^{j+1}\le n$, compare $N$ and $n$ modulo $p^{j+1}$. Otherwise $p^{j+1}>n$ gives $v_p(n)=j$ directly. Put $v=n/U$. We have positive rows of weight $n$,
$$
 A_0=(a,b,c+1),\qquad C_0=v(f,g,h).
$$
The parts of $A_0$ are divisible by $U$, those of $C_0$ by $v$, and $v\ge3$.

For $p\ne q$, the parts of $A$ reduce to those of $A_0$ at all layers $p^j\le n$. If $p\nmid U$, this holds for $C,C_0$ as well, by dividing $N\equiv n$ by the unit $U$. A carry in a row at weight $n$ has a witness at one of these layers: the sum of its three least residues is greater than $n\bmod p^j$. It therefore persists at the large weight for the rows just compared.

It follows that $A_0$ avoids the primes of $U$ and $C_0$ avoids the primes of $v$ with the possible exception of $q$. At every prime $p\mid n-1$, one of the small rows avoids $p$. For $p\ne q$ this transfers from $N-1$, since $p\nmid U$. For $p=q\mid d$, the slopes $a,b,c$ add without a base-$q$ carry: the large row has precisely their digits above its bottom $e$ places. Their units digits must then all vanish because $q\mid d$. Adding one to $c$ leaves their addition carry-free. Finally, if the assigned row $C_0$ carries at $q\mid v$, **decreasing the distinguished part of the opposite row $A_0$ by one gives a positive carry-free triple at $q$**, again by the digits above the bottom $e$ places in $A$. This last test is essential in the computation.

## Finite test at the small weight

Assign the full powers of $n-1$ to the two avoiding rows, making products $s,t$. For the two weight factors rename the scales $u,v$ so that the row scaled by $u$ takes $s\ge v$. Such an orientation always exists since $uv=n$ and $st=n-1$. The analysis at the [endpoint](endpoint_reduction.md#splitting-the-two-consecutive-integers) applies even with the possible missing test at one prime of $n$: divisibility of **all** parts by their scale still holds. Each power assigned from $s$ determines two zero residues and one exceptional residue in the normalized parts modulo that power. The positive normalized parts are all below $v\le s$, so each of the three positions must be exceptional for some prime. In particular $s$ has at least three prime factors. This proves the assertion for $\omega(d)\le2$. If $n$ is a prime power or twice an odd prime power, it has no unitary factorization with both factors at least 3, proving the other elementary cases.

For $3\le n\le10^9$, the C++ verifier enumerates all size-feasible factorizations as in the [original range computation](full_range_1b.md). It checks the larger-modulus row by the Chinese remainder theorem and digit addition, allowing it to miss at most one assigned prime of the weight. Among 6,151,200,914 size-feasible assignments, it emits 236 rows. In order, the counts in blocks of $10^8$ are

~~~~
26 29 16 23 17 21 27 31 19 27
~~~~

The Python verifier refactors the assigned products by trial division and validates every emitted row by factorial valuations. It enumerates **all** positive unordered triples of the complementary normalized sum, at most 179. Across 10,554 companion rows, ten pairs pass the carry tests when one assigned weight prime is waived. In each of the ten, all three ways of decreasing a part of the opposite row fail the final test above. Thus no pair survives.

All ten waived primes equal 2. The sorted residues of the opposite row modulo eight are $(2,5,5)$ five times, and $(2,7,7),(3,6,7),(4,5,7),(4,6,6),(5,5,6)$ once each. Each sum exceeds eight, so after decreasing any part by one the residue sum is at least eight. In the original large row this forces a carry at the layer $8Q$ for every binary scale $Q=2^e$. The [Lean proof](WeightTransfer.lean) verifies the residue inequalities and the implication for every scale $Q\ge2$ (a power of two makes $8Q$ a binary layer).

## A density-one set of unbounded cofactors

The count for large $d$ uses the part congruences on both rows and the binary digits of just one of them. If $q\ne2$ or $2\nmid n$, the small rows include a binary carry-free row. Assign 2 to such a row in the appropriate factorization (if $2\mid n$, the assignment already supplies it). The proof of the [bound of exponent 0.8689](two_predecessor_congruences.md#the-unrestricted-three-row-count) uses only these binary digits, the row divisibilities by the two scales and the three residue positions at every predecessor power. All those hypotheses have just been proved, even when the test at an odd weight prime is lost. Consequently it counts the union over these $q$.

For $q=2\mid n$, use $A_0$ and let $S$ be its assigned predecessor product. Now $U,S$ are both odd. Its parts $x_i$ obey $U\mid x_i$ and $S\mid x_i(x_i-1)$, and each $x_i\ge3$. The parts $y_i=(a,b,c)$ partition the binary ones of $d=n-1$. Each $y_i$ determines $x_i$ up to adding zero or one, so there are still $X^{o(1)}$ possible pairs $(U,S)$ given any one $y_i$ when $d\le X$.

To spell out the binary count, put $K=US$, and take a dyadic block of $d$ of bit length $l$ with $X=2^l$. For $K\le X^k$, the [two-congruence reconstruction](two_predecessor_congruences.md#recovering-the-row-from-a-product) on the complementary row bounds the number by $X^{k+o(1)}$ (the boundary weight $n=X$ can be counted separately). Otherwise a specified pair $(U,S)$ fixes $d$ in one class modulo the odd $K$. Either the first or the last $h=l-\lfloor kl\rfloor$ bits of $d$ determine $d$ in that class. The first bound uses an interval of length $2^{l-h}\le K$ and the last uses modulus $2^h K\ge X$. Choose a canonical witness for each $d$. The [binary entropy bound, equations (3) and (4)](common_factor_count.md#information-that-every-individual-part-supplies), now applies to $y_1,y_2,y_3$ summing to $d$, with no fixed low places to strip: each part determines the modulus factors up to $X^{o(1)}$ choices, and their bits partition those of $d$. Taking $k=0.86884$ and the checked parameters in the exponent-0.8689 bound gives $O(X^{0.8689})$ after summing blocks. This counts all possible exceptions in the remaining case.

## Reproduction

With a C++17 compiler supporting 128-bit integers and Python 3.9 or later, run from the repository root.

~~~~sh
mkdir -p .work
c++ -std=c++17 -O3 check_endpoint_range.cpp -o .work/transfer_range
.work/transfer_range --transfer 3 1000000000 > .work/transfer_rows.txt
shasum -a 256 .work/transfer_rows.txt
python3 check_weight_transfer.py < .work/transfer_rows.txt
python3 check_weight_transfer.py --controls
lean WeightTransfer.lean
~~~~

The SHA-256 of the sorted-by-weight enumeration output in its original within-weight order is `bc0da97d6580dba9293fc9dd752e236bb1e9dbb01d8a642bd41b83a08b3860af`. Independent Python assignment generation covers the first 20,000 weights and recovers the first emitted row at $n=8{,}284{,}703$, where the waived prime is 19. The controls also compare the carry and one-part decrement tests with exact binomial products through weight 74.
