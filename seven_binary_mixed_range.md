# A binary base and an odd predecessor base, both with at most seven units

**Finite theorem.** Let $3\le N<2^{1,000,000}$. If $2$ and an odd prime $p$ divide $N-1$ and $s_2(N),s_p(N)\le7$, Wasserman's three-row conjecture holds at $N$.

**Upper cutoff.** The same assertion holds for $N\ge2^{10^{180}+1}$. Together with the [all-weight result for two odd primes](seven_mixed_primes.md), this leaves only a finite intermediate range for the assertion with *any* two distinct predecessor primes of digit sum at most seven. The intermediate range has not been exhausted.

The [five-unit theorem](six_digit_gaps.md) deals with either sum at most five. Since $N$ is odd, the remaining odd-base sum must be seven. By the [large-base theorem](seven_large_bases.md), only $3\le p<257$ needs to be considered. Put $n=N-1$. It has five or six binary ones, all at positive positions, and is the sum of six positive powers of $p$, with multiplicities less than $p$.

## The lower range

First enumerate $n\le2^{128}-2$. Choose the six positive $p$-powers in nonincreasing order. At a partial sum $z$ with $k$ places left, the final sum is at least $z+kp$. For any positive integer $w$, the least integer $w'\ge w$ with at most six binary ones is $w$ if it already has at most six. Otherwise keep the six highest ones, erase the lower ones and add the value of the sixth one. Indeed all integers from $w$ up to but excluding this next multiple of that bit still have the same six leading ones and a nonzero tail. The multiple has at most six ones after carrying. Use this operation on $z+kp$. The next summand must be at least $(w'-z)/k$, which permits a binary search in the list of powers. The verifier enforces the multiplicity bound and tests completed sums exactly.

A second enumeration starts from all five- and six-element subsets of the positive bit positions below 34 and tests ordinary digits in each eligible odd prime base. It recovers the entire same list, including the maximum $N=10,737,942,553$ at $p=3$. There are 1,987 weight/base pairs. Each is certified using the [seven-unit row bound](seven_digit_units.md#the-inner-units). Enumerating all 299,806 positive $p$-avoiding rows leaves 380 row gcds $U_0\ge7$. For each, the odd part of the product over all part pairs of $\gcd(x_i,x_j,n/p^{v_p(n)})$ is an upper bound $S_0$ on the assigned odd predecessor factor. Every row fails $27(n/p^{v_p(n)})^2\le U_0^3 S_0^2$. [Exact lower-range verifier](check_mixed_sparse_range.py)

For larger weights the integer prefix sieve considers only the leading power $p^e$. For a leading power smaller than $2^{120}$ the entire weight is already covered, since $n\le6p^e$. For every other possible leading power below $2^{1,000,000}$:

1. Enclose $p^e$ between $l2^s$ and $u2^s$. Starting with $l=u=1,s=0$, multiply both endpoints by $p$, then shift the endpoints right as needed to retain 116 significant bits, rounding the lower endpoint down and the upper endpoint up. Accumulate the shifts in $s$. All products and rounding additions fit in 128 bits since $p<256$. Terminate only when the **lower** endpoint reaches the specified ceiling.
2. Round outward again so that about 100 leading bits are retained. At this fixed binary scale divide the two endpoints successively by $p$ with outward rounding. This list encloses every available smaller power; the final interval $[0,1]$ includes all omitted smaller powers.
3. If a completion exists, the floor of its scaled sum has at most six ones: it is a high binary prefix of $n$. Search over the five remaining scaled summands in nonincreasing order, allowing even excess multiplicities. At a partial interval $[l',u']$ with $k$ places remaining, let $w'$ be the first sparse integer at least $l'$. If $w'\le u'+k$, retain the leading power as a possible hit. Otherwise only summands whose upper endpoint is at least $(w'-u')/k$ can be the next one. Exhausting these choices rejects the leading power.

All 9,391,604 tested leading powers were rejected, using the 53 eligible odd primes. Positive controls normalize known small solutions to the large mantissa range and require the search to retain them. The verifier also compares the sparse-successor operation with direct counting. An additional full-integer enumeration through 32,768 bits, using the lower-range algorithm with the higher ceiling, gave the same 1,987 pairs. [Integer prefix sieve](check_mixed_sparse_prefix.cpp)

Reproduce the two principal checks in the repository root:

~~~~sh
python3 check_mixed_sparse_range.py
mkdir -p .work
c++ -std=c++17 -O3 check_mixed_sparse_prefix.cpp -o .work/mixed_sparse_prefix
.work/mixed_sparse_prefix 1000000
~~~~

For the longer full-integer comparison, use `python3 check_mixed_sparse_range.py --bits 32768`.

## An absolute upper cutoff from the two bases

The joint boundedness of digit sums in fixed multiplicatively independent bases is the effective theorem of Stewart. Here the reduction to finitely many primes gives an absolute bound. We give its size explicitly using the same logarithm estimate as in the [two-odd-prime proof](seven_mixed_primes.md). [Stewart, On the representation of an integer in two different bases](https://uwaterloo.ca/pure-mathematics/sites/default/files/uploads/documents/j-reine-ange-math-1980.pdf)

Write the binary exponents of $n$ in decreasing order $E=e_1>e_2>\cdots$, and the occupied $p$-positions as $D=d_1>d_2>\cdots$ with their actual coefficients. Each side has at most six units. Its prefix down to $e_i$ or $d_j$ has the form $2^{e_i}R_i$ or $p^{d_j}T_j$, respectively. Define
$$G=\max((E-e_i)\log 2,(D-d_j)\log p).$$
Then $R_i,T_j\le6e^G$. If the two prefixes coincide, $2^{e_i}\mid T_j$ since $p$ is odd, so $E\log 2\le2G+\log6$.

Otherwise put $\Gamma=2^{e_i}p^{-d_j}R_i/T_j-1$. Matveev's theorem for positive rationals gives
$$\log|\Gamma|>-10^{13}(1+\log(E+1))(G+2).$$
The heights of $2,p$ are less than 1 and 6, and the height of the reduced rational $R_i/T_j$ is at most $G+\log6$. The exponents are bounded by $E+1$ because $p^D\le n<2^{E+1}$. The three-logarithm coefficient $1.4\cdot30^6\cdot3^{4.5}\cdot6$ is below $10^{13}$. If the rational factor equals 1, use the two-logarithm bound. [Matveev's rational bound, Theorem 2.1](https://archive.mpim-bonn.mpg.de/id/eprint/5025/1/mpim-preprint_2023-34.pdf)

Let $g$ be the smaller of the two gaps from the leading exponent to the next unused position, in the same logarithmic units, ignoring an exhausted side. Each unused tail has size at most $6n e^{-g}$. The denominator prefix is at least $p^D\ge n/6$. Hence $|\Gamma|\le36e^{-g}$. Extend the prefix on the side giving $g$. For the new maximum gap $G'$ these bounds imply
$$G'+2 \le C(G+2),\qquad C=10^{14}(1+\log(E+1)).$$
Starting with just the two leading positions gives $G=0$. After at most ten extensions the prefixes agree, since both contain their respective entire expansions by then. Thus $E\log2\le4C^{10}$ and $E<6C^{10}$. For $E=10^{180}$ the right hand side is already smaller than $E$, and $(1+\log(E+1))^{10}/E$ decreases thereafter by differentiation. In particular $E<10^{180}$, proving the stated upper cutoff.
