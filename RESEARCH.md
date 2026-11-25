# Prime certificates for Wasserman’s trinomial conjecture

Research checkpoint, 25 November 2026. A certificate for a weight $N$ is a set of primes from which any three positive trinomial coefficients of weight $N$ have a common divisor. The prime may depend on the three coefficients. Repeated triples are allowed.

The [later finite theorem](full_range_1m.md) extends coverage of the full three-row assertion to one million. The certificate-size results below concern weights through 5,000.

## An independently checked finite theorem

**Every weight from 3 through 5,000 has a certificate of size at most five.** The minimum is one for 1,788 weights, four for 3,208 weights, and five for exactly two weights, 699 and 717. The table gives one minimum-size certificate per weight.

The Python generator uses Kummer’s criterion. Write $F_p$ for the set of unordered positive triples whose coefficient is coprime to $p$. To test a proposed set $S$, assign its primes to three groups. A group is feasible precisely when some row belongs to every $F_p$ indexed by it. Enumerating the digit allocations in one base and testing the remaining bases decides feasibility. Empty groups may be assigned any positive triple. An assignment of all primes to feasible groups is equivalent to three rows collectively avoiding every prime of $S$. Given such an assignment, the generator adds a prime dividing all three coefficients. When there is no feasible assignment, the current set certifies the weight.

When the search returns more than four primes, the generator also searches all possible four-prime certificates using the disjointness criterion proved below.

The standalone C++ verifier takes the prime sets as untrusted input. It checks primality and consecutive coverage, loops over every $1\le a\le b\le c$ with $a+b+c=N$, and computes prime valuations from Legendre’s factorial formula. For each observed mask of prime divisors it checks all three-way intersections, including repetitions. A separate Python verifier proves the size lower bounds. For 699 and 717 it evaluates a fixed sample of coefficients by exact binomial products and excludes every four-prime set on that sample.

Reproduce in this directory with Python 3.9+ and a C++17 compiler.

```sh
python3 check_wasserman_785.py
python3 prime_certificates.py 3 5000 > regenerated.txt
cmp regenerated.txt certificates_5000.txt
python3 check_minimality.py
c++ -std=c++17 -O3 verify_certificates.cpp -o verify_certificates
./verify_certificates 3 5000 < certificates_5000.txt
echo '785 2 3 5 7' | ./verify_certificates 785 785
```

The interval check counted **3,473,263,472** unordered positive triples. The SHA-256 of the supplied table is `f02c33b5589456b82fefa79bc548606940cdbb4cab2c7a64632bfe9df9339132`.

At 785 the supplied exact binomial-product script and the independent verifier give the following frequencies.

| gcd with 210 | 30 | 42 | 70 | 105 | 210 |
|---|---:|---:|---:|---:|---:|
| Unordered positive triples | 11 | 16 | 11 | 6 | 51,308 |

Every coefficient misses at most one of 2, 3, 5, 7. Three rows can therefore miss at most three of these primes between them.

## The first weight requiring five primes is 699

All weights below 699 have a witness of size at most four in the verified table. The set $\{2,41,233,349,691\}$ certifies 699. Here is the matching lower bound.

A prime divides every positive trinomial of a given weight precisely when the base-$p$ digit sum of the weight is at most two. In one direction, a carry-free sum of three positive integers uses at least three digit units. In the other direction, distribute the units of the digits among three nonempty boxes. The distribution produces a carry-free triple. For odd $N$, a sum of two powers of an odd prime is even. Consequently an odd $N$ can have a universal prime only if it is a prime power or has the form $2^a+1$. Since $699=3\cdot233$ and $698=2\cdot349$, no prime is universal at 699.

With no universal prime, a four-prime set certifies a weight exactly when its four sets $F_p$ are pairwise disjoint. A row in two of them and one row from each of the remaining two would defeat the certificate. Conversely, three rows cannot meet four disjoint sets. The same argument rules out any set of at most three individually avoidable primes.

At 699, the row $(1,1,697)$ has coefficient $699\cdot698$. Thus at least three primes of any four-prime certificate must belong to $T=\{2,3,233,349\}$. For each of $p=2,3$, the rows below have coefficients coprime to $p$.

| $p$ | Rows | Exact gcd of the row coefficients |
|---|---|---:|
| 2 | (1,2,696), (1,42,656), (1,18,680), (3,48,648) | $81317=233\cdot349$ |
| 3 | (3,3,693), (3,279,417), (3,174,522), (84,246,369), (3,36,660), (18,324,357), (3,48,648) | $28216999=233\cdot347\cdot349$ |

If 2 belonged to a certificate of size four, all three other primes would divide 81,317 by pairwise disjointness, which is impossible. Without 2, the certificate must contain 3, 233 and 349. The second gcd forces 347 as its fourth member. But $(1,349,349)$ is carry-free in bases 347 and 349. In the first base its units sum to 5 and its next digits to 2. This proves the lower bound.

To reproduce either gcd with exact integers, set `rows` to the corresponding list and run

```python
from math import comb, gcd
values = [comb(699, a) * comb(699-a, b) for a, b, c in rows]
assert all(a+b+c == 699 and a > 0 and b > 0 and c > 0 for a, b, c in rows)
print(gcd(*values))
```

At 717 the certifying set is $\{2,13,179,239,709\}$. The independent exact check uses $717\cdot716$ and the same disjointness argument. Its witness rows are included in the verifier.

## Literature checked

Bergman’s Proposition 15 gives the strict bound $N<785$. Proposition 10 excludes a row $(N-i-j,i,j)$ with $3\le i+j\le10$ whenever three rows collectively avoid every prime dividing $N(N-1)(N-2)$. The addenda work through 78 and give additional restrictions in Section 6. The author’s publication list, searches for the conjecture by name, and forward citation searches in OpenAlex and Semantic Scholar did not reveal a later published improvement to the consecutive-weight bound. The finite theorem above is a project computation. Its priority is undetermined. [Published paper](https://doi.org/10.1017/S0004972710001723) · [arXiv](https://arxiv.org/abs/0806.0607) · [Addenda](https://math.berkeley.edu/~gbergman/papers/unpub/from_nomial.pdf) · [Author’s publication list](https://math.berkeley.edu/~gbergman/papers/)

Casacuberta studies a related prime-covering problem, including multinomials in Section 6. In that problem each individual coefficient must be divisible by some member of a fixed set. Bergman revisits the orbit question as Q11 in his [later notes](https://arxiv.org/abs/1904.04298). Rukavicka obtains divisibility under hypotheses on the multiplicities of an ordered rooted tree. Neither supplies the simultaneous common prime asserted here. [Casacuberta](https://arxiv.org/abs/1906.07652) · [Rukavicka](https://doi.org/10.1007/s00013-015-0770-5)

## An explicit infinite family from Proposition 10

For $e\ge4$ and $0\le t\le10$, let $N=2^e+t$. In a putative counterexample one row must be carry-free in base 2. Its largest part receives the $2^e$ digit, so its two smaller positive parts sum to at most $t$. Their sum cannot be 2, since $1+1$ carries in binary. Therefore the sum lies from 3 through 10, which Proposition 10 excludes. This deduction from Bergman gives all eleven offsets, uniformly in $e$. For a subfamily without any universal prime, take $e=20k+5$ and $t=3$ for $k\ge0$. Then $N\equiv10\pmod{25}$, $N>5$, and $N-1=2^e+2$ is not a power of two.

A separate argument in [Infinite families with bounded prime certificates](bounded_certificates.md) gives a uniform bound of seven primes for infinitely many odd composite weights without a universal prime.
