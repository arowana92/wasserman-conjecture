# Every nonendpoint part of a counterexample exceeds one million

In any counterexample the [endpoint theorem](all_weights_endpoint_reduction.md) requires one row to be a permutation of $(1,1,N-2)$. **Every part of either of the other two rows exceeds $1{,}000{,}000$, at every weight $N$.** The two finite enumerations below have no input cutoff on the weight. The largest weight reconstructed before the carry tests was $191{,}999{,}809{,}000{,}000$.

**Intervals above sums of two powers.** At every prime $p$ and every weight $N$ with $p\mid N(N-1)$, delete the two largest units of its base-$p$ expansion, or all units if there are fewer than two. Write $r$ for the sum that remains. If $r\le1{,}000{,}000$ or $r^5\le27(N-1)^2$, the full conjecture holds. In particular for every $i,j\ge0$ the entire integer interval
$$[2^i+2^j,\ 2^i+2^j+1{,}000{,}000]$$
satisfies the conjecture at weights at least three.

Indeed, a nonendpoint row avoiding $p$ must allocate the digit units without a carry. At least one of its three positive parts receives neither deleted unit and is at most $r$. Apply the new lower bound or the [individual-part inequality](endpoint_factors.md). For the interval claim, carrying in an expression as a sum of powers replaces several units by their sum. It never decreases the sum of the two largest units, padding with zeros when necessary. Thus at $N=2^i+2^j+d$ the remainder is at most $d$, even for repeated powers or overlapping low digits.

## Enumeration without a cutoff on the weight

Let $a$ be any part of row $A$. In a counterexample assign the full prime powers of $N$ and $N-1$ as $UV=N$ and $ST=N-1$, with $U,S$ assigned to $A$ and $V,T$ to row $C$. The coefficients avoid every prime assigned to their row. Every part of $A$ is a multiple of $U$, and the three positive parts of $C/V$ sum to $U$. Thus $U\ge3$ and $a=Ux$ for some positive $x$.

At each full power of $S$, the actual residues of $A$ are $0,0,1$ in some order. Let $F$ be the product of those at which $a$ is zero and $B$ the product at which it is one. Then
$$F\mid x,\quad B\mid a-1,\quad S=FB,\quad k=(a-1)/B.$$
The other two parts of $A$ are positive multiples of $UB$. The three parts of $C/V$ have product divisible by $T^2$, and two are positive multiples of every full power of $T$. Hence
$$27T^2\le U^3,\qquad 2w+1\le U\quad(w\parallel T\text{ a full prime power}).$$
Finally $FT\equiv k\pmod U$, by $ST+1=UV$ and $kB=a-1$. In particular a fixed $a$ gives a finite list of weights and assignments.

The first verifier loops over $a,U,F,B$ and the single progression for $T$. It requires coprimality of the factors on each side of $UV$ and $ST$. It checks that $N-a\ge2UB$ and that $\binom Na$ is coprime to $US$. The latter is necessary since $\binom Na$ divides the full coefficient of $A$.

Exactly one of $T>U$ and $S>V$ holds. This follows from $ST=UV-1$, coprimality of $N$ and $N-1$, and $U,V\ge3$. If $T>U$, each part of $C/V$ is a positive integer below $T$. At $w\parallel T$ its residue is $0$ or $U$ since $VU\equiv1\pmod w$. Assign each power to its unique nonzero position. All three positions must be used. The Chinese remainder theorem fixes each part from the assignment. Check the sum and test that $C$ avoids all primes of $T$.

If $S>V$, write the other two parts of $A$ as $UBh,UB(l-h)$, where $l=(V-x)/B<F$. At each full power of $F$, one of $h,l-h$ vanishes. Enumerate the residues $h\equiv0$ or $l$ at each such power. Retain $1\le h\le l/2$ and check the row in every base dividing $US$.

The independent verifier starts from $U,x$ and the divisors $S$ of $x(a-1)$. It enumerates $V\equiv U^{-1}\pmod S$ directly using $V\ge x+2\gcd(S,a-1)$ and $27((UV-1)/S)^2\le U^3$. For $S>V$ it tests the possible splits by scanning the multiples of $UB$. For $T>U$ it computes admissible single entries by the unitary divisors of $T$ and looks for three that sum to $U$. It takes successive factorial floors for each row test. Both code paths count the assignment for each specified value of $a$ and allow repetitions in the row.

| Stage | Both implementations |
|---|---:|
| Assignments after the size bound on $T$, the lower bound on $N-a$ and unitarity | 2,453,363,033 |
| Cases $S>V$ after the binomial test and individual-power bounds | 152,853 |
| Cases $T>U$ after those tests, with at least three prime bases in $T$ | 35,715,732 |
| Rows $A$ passing the required bases in the first case | 0 |
| Rows $C$ passing the required bases in the second case | 0 |

As a positive control both row routines reconstruct $(45{,}491{,}401,86{,}617{,}810,90{,}679{,}143)$ at $N=222{,}788{,}354$. They recover it from each of its three parts in the first routine and from the complementary scale alone in the other. Its zero valuations at $29,127,241,283$ are checked separately in the [three-factor study](endpoint_factors.md). The control tests the row routines in isolation.

~~~sh
mkdir -p .work
c++ -std=c++17 -O3 check_fixed_parts.cpp -o .work/check_fixed_parts
c++ -std=c++17 -O3 check_fixed_parts_reverse.cpp -o .work/check_fixed_parts_reverse
.work/check_fixed_parts 3 1000000 0
.work/check_fixed_parts_reverse 3 1000000
.work/check_fixed_parts --example
.work/check_fixed_parts_reverse --example
~~~
