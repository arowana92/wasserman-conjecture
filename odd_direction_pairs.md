# All common weights for primitive row sums through 160,000

**The three-row conjecture holds at every weight if at least two of the three rows, each divided by its own coordinate gcd, have sums at most 160,000.** There is no bound on the third primitive sum or on the original weight.

For the underlying two-row result, take any two primitive positive triples with sums $u,v\le160000$. Scale them to any common weight $N$. Their coefficients share a prime dividing $N(N-1)$. The [unique-scale theorem](primitive_row_scales.md) settles every weight except a possible single weight per pair, where all lowest-layer avoidance tests pass. On the entire list of those possible weights, **six primes suffice and six are necessary**. One smallest certificate is
$$
                           \{2,3,5,7,11,23\}.
$$
The prime used may depend on the pair. At scales already settled by the lowest layers, the prime can also depend on the scale. For the three-row statement, a failure must [contain the endpoint](all_weights_endpoint_reduction.md). At $N>160000$ its primitive sum is $N$, so the two bounded directions must be the other rows, excluded by the two-row result. The [verified weight range](full_range_1b.md) covers the remaining $N$.

## Enumerating the pairs that need a prime other than 2

Write $N=guv$, $\gcd(u,v)=1$, and assign the full powers of $N-1=ST$ to the two primitive rows $x,y$ passing their lowest layers. At least one row, say $x$, has $S\ge u$. This follows from $ST\ge uv-1$. Group the powers of $S$ by their exceptional coordinate. No group is empty since each positive entry of $x$ is less than $u\le S$. The three group products are pairwise coprime, so after permutation they are $2\le d<e<f$. The three entries of $x$ are positive multiples of $ef,df,de$ respectively, with sum at most $\min(def,X)$ where $X$ is the input bound. The program enumerates these integers directly and checks primitivity.

For each $x$, the multiplier is $m=N/u=(ST+1)/u$, hence $T\equiv-S^{-1}\bmod u$. We need only test
$$
                  1\le T\le\sqrt{X^3/27},\qquad m\le X^2/27+1.
$$
Indeed $27S^2\le u^3$ and $27T^2\le v^3$ by the product of the primitive entries. A parity survivor must have an odd coefficient on some row. If $T\ge v$, both rows meet the first-row condition, so it suffices to generate that assignment with an odd first coefficient. In particular, if two entries of $x$ have the same 2-adic valuation, their lowest occupied bits coincide at all scales and one need only enumerate $T<X$. In the other cases the actual disjoint-bit test on $mx$ decides its parity.

Discard $\gcd(S,T)>1$. Try each $v\mid m$ with $3\le v\le X$, $\gcd(u,v)=1$, and $27T^2\le v^3$. At a power $w$ assigned to a coordinate of $y$, that coordinate must be $v\bmod w$ and the other two must vanish. Enumerate the unordered partitions of the full powers of $T$ into up to three groups. The CRT specifies all three residues modulo $T$. Enumerate the positive representatives with sum $v$, require primitivity, and retain the pair if either scaled row has odd coefficient. If the companion has to be odd then its primitive sum $v$ is odd, since its primitive entries have a unique least binary valuation, equal to zero. The code deduplicates assignments and ordering. These steps give the [complete list of 272 pairs](odd_direction_pairs_160000.txt).

| Primitive sum bound | Pairs with an odd coefficient | Distinct candidate weights with such a pair | Largest of those weights |
|---:|---:|---:|---:|
| 6,000 | 5 | 3 | 13,073,521 |
| 12,000 | 17 | 13 | 63,422,464 |
| 30,000 | 39 | 33 | 457,699,725 |
| 100,000 | 183 | 125 | 16,819,136,685 |
| 160,000 | 272 | 195 | 67,104,656,541 |

For every record, Python independently checks primitive-scale recovery, the actual prime-power residues of both rows, and the parity. It factors $N(N-1)$ and computes the common endpoint gcd from Legendre's factorial floors, checking each valuation by the base-$p$ digit-sum formula. The smallest gcd on the 272-record list is 675,213. A separate [enumeration with no parity pruning](counting_primitive_directions.md) returns exactly the same records through sum 6,000.

## Sharp changes in the size of the certificate

The table concerns a single set of primes covering **every pair that passes the lowest layers**, evaluated at each pair's candidate weight. The first nonempty list occurs at sum bound 149, with directions $(4,17,68),(39,45,65)$ and weight 13,261.

| First sum bound | Minimum certificate size | A certificate until the next change, or through 160,000 |
|---:|---:|---|
| 149 | 1 | $2$ |
| 2,769 | 2 | $2,3$ through 7,150 and $2,5$ through 8,191 |
| 8,192 | 3 | $2,3,5$ |
| 35,671 | 4 | $2,3,5,7$ |
| 96,773 | 5 | $2,3,5,7,11$ |
| 100,547 | 6 | $2,3,5,7,11,23$ |

At the last change the weight is $4,019,768,513$, with primitive directions $(643,41152,58752)$ and $(2852,3151,33976)$. Their common gcd with the endpoint is exactly $62,808,883=23\cdot31\cdot137\cdot643$.

For a certificate containing 2, its remaining primes must meet the common-prime set of each parity survivor. The verifier computes the exact minimum by selecting an unhit set and branching over each of its primes. Six explicit first-layer witnesses have pairwise coprime odd parts of their endpoint coefficients and primitive sums below 2,769. A certificate omitting 2 therefore requires at least six primes at every transition from 2,769 onward.

There is a short factorization-free check of the sharp six-prime lower bound, independent of the large enumeration. [Lean](DirectionPrimeCovers.lean) checks seven first-layer witnesses with an odd coefficient and primitive sums at most 100,547. Their seven odd endpoint parts require at least five divisors greater than one to hit them. It also checks the six witnesses for a certificate omitting 2. For the covering proof, put each integer into a bin, replacing the bin value by its gcd with that integer. A covering list of divisors supplies a placement in the same number of bins with no gcd equal to one. Lean proves this implication for every list, then exhausts the placements for the explicit integers by kernel reduction. For a certificate containing 2, five further primes are necessary because 2 misses all seven odd-coefficient witnesses. This Lean file verifies the lower bound and the witness data. The upper bound and enumeration use the C++ and Python programs.

## Reproduction

From the repository root, with a C++17 compiler and Python 3.9 or later, run

~~~sh
mkdir -p .work
c++ -O3 -std=c++17 check_odd_direction_pairs.cpp -o .work/odd_directions
python3 check_odd_direction_pairs.py .work/odd_directions --limit 160000 --jobs 8 --compare 6000 --output .work/odd_pairs.txt
cmp .work/odd_pairs.txt odd_direction_pairs_160000.txt
lean DirectionPrimeCovers.lean
~~~

The final Lean check was run with Lean 4.36.0. The sorted record file has SHA-256 digest
~~~text
c0cb511c2cd46efffe85a546943bcde3ffa719f3e3b01741ad55a03f8f691eab
~~~
