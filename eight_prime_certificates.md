# Weights requiring eight primes

For fixed $N$, a row is a positive triple of sum $N$. A prime list certifies triples of rows if for every three rows some prime on the list divides all three trinomial coefficients. A list certifies pairs if the corresponding statement holds for every two rows. A row *avoids* a prime if its coefficient is a unit at that prime. Write $D$ for the prime divisors of $N(N-1)$, the coefficient of the endpoint row $(1,1,N-2)$.

| Weight $N$ | Minimum for pairs using $D$ | Minimum for triples using any primes | Minimum for triples using $D$ |
| ---: | ---: | ---: | ---: |
| 12,189,629 | 5 | 6 | 7 |
| 7,064,355,413,422,079 | 7 | 8 | 10 |
| 16,795,978,692,415,487 | 7 | 8 | 10 |

At the second weight one minimum list consists of $179,181,271,283,1129,1627,7487$ and $N-22$. At the third use $193,409,691,2447,5147,6469,6779$ and $N-30$. Thus any proposed uniform bound for the number of primes in a certificate for triples of rows must be at least eight. These are pointwise results about all rows at the stated weights.

## From pairs to triples

Suppose $D$ certifies pairs and let $k>1$ be the size of its smallest subset certifying pairs. Equivalently, for every subset of fewer than $k$ endpoint primes there are two rows which between them avoid the whole subset. Each individual prime can be avoided by a row.

**Lower bounds.** Every certificate for triples has at least $k+1$ primes. If a proposed list of size at most $k$ contains an exterior prime, use the endpoint row to avoid every exterior prime. Two more rows can avoid the remaining primes. For a nonempty list contained in $D$, reserve one prime for one row and avoid the rest with two rows.

A list $S\subseteq D$ certifying triples has size at least $k+\lceil(k-1)/2\rceil$. For any row, the primes of $S$ dividing its coefficient must certify every pair, so there are at least $k$ of them. In particular $|S|\ge k$. Pick $k-1$ members of $S$ and two rows avoiding them between them. At least one row avoids $\lceil(k-1)/2\rceil$ of these primes, giving the bound.

**Equality for unrestricted lists.** For $k\ge4$, a list of size $k+1$ certifies triples precisely when it has the form $C\cup\{r\}$, where $C\subseteq D$ certifies pairs, $|C|=k$, $r\notin D$, and every row avoiding $r$ is divisible at every prime in $C$. Indeed, the previous bound excludes lists contained in $D$. The endpoint row forces their intersection with $D$ to certify pairs. This leaves exactly one exterior prime. If a row avoiding $r$ also avoids a member of $C$, two further rows can avoid the remaining $k-1$ primes. Conversely, given three rows, either $r$ divides all their coefficients or one row avoids it. Choose a common prime of $C$ for the other two rows. It divides the first coefficient as well.

**Which exterior primes can occur?** Form the multiset of base-$p$ digit units of $N$ by including $p^i$ with multiplicity equal to digit $i$. For $p\in D$ let $h_p$ be the sum of the two least units. There are at least three units because $k>1$ and some positive row avoids $p$. By Kummer’s criterion the parts of a row avoiding $p$ partition this multiset. The smallest attainable sum of its two smallest parts is $h_p$, attained by putting the two least units in separate parts and all remaining units in the third.

For a prime $N/2<r\le N$ exterior to $D$, a row avoids $r$ precisely when its two smallest parts sum to at most $N-r$. Thus the incompatibility required for $p\in C$ is exactly $h_p>N-r$.

There is a finite restriction covering smaller exterior primes as well. Let $G$ be the $k$th largest value among the $h_p$, including ties. Every exterior prime in a size-$k+1$ certificate divides some $N-j$ with $0\le j<G$. Indeed, $C$ contains a $p$ with $h_p\le G$. Use the $p$-avoiding row built from its two least digit units. Its coefficient divides $N(N-1)\cdots(N-h_p+1)$ and must be divisible by $r$.

## How the numbers are checked

At $N=12{,}189{,}629$, the endpoint factorizations are $N=17\cdot53\cdot83\cdot163$ and $N-1=2^2\cdot11\cdot29\cdot41\cdot233$. The 34 supplied rows suffice to avoid every four endpoint primes in pairs. Among $17,29,41,53,83,163,233$ no row avoids three, so any five certify pairs and all seven certify triples. The nine rows avoiding the prime $N-6$ have coefficients divisible by all seven endpoint primes. For example, $17,53,83,163,233,N-6$ certify triples. [Witnesses](six_prime_witnesses.txt) · [Verifier](check_six_prime_minimum.py)

At the two larger weights the factorizations are squarefree.

| $N$ | Factors of $N$ | Factors of $N-1$ |
| --- | --- | --- |
| 7,064,355,413,422,079 | 67, 101, 103, 127, 181, 271, 1627 | 2, 73, 113, 179, 283, 1129, 7487 |
| 16,795,978,692,415,487 | 101, 127, 193, 409, 2447, 6779 | 2, 53, 71, 97, 691, 5147, 6469 |

For these weights, respectively, 189 and 142 supplied rows give two avoiding rows for each of the 3,003 and 1,716 six-subsets of $D$. Their valuations are checked by Legendre’s formula and independently by digit addition. In each of the displayed seven-prime lists no four primes are jointly avoidable. Two rows can therefore avoid at most six members, proving $k=7$. The near-prime test above proves the eight-prime upper bounds. Both near primes have full Lucas primality certificates. [Witnesses](eight_prime_witnesses.txt) · [Verifier](check_eight_prime_minima.py)

A ten-subset of $D$ certifies triples precisely when no four of its members are jointly avoidable. If four can be avoided by one row, use the six-subset witnesses for the complement. Conversely, three rows avoiding at most three each leave at least one prime. At the first large weight the only jointly avoidable four-subset is $\{2,179,1129,7487\}$. There are $\binom{14}{10}-\binom{10}{6}=791$ minimum endpoint-only lists. At the second large weight, 20 of the 715 four-subsets are jointly avoidable, leaving 33 minimum ten-subsets. One is $\{53,71,127,193,409,691,2447,5147,6469,6779\}$.

At the first large weight $G=254$. Each exterior prime $r>N/2$ in a minimum list therefore has $2\le N-r<254$. Every such prime in the interval has $N-r\ge3=h_2$. Consequently any seven eligible primes exclude 2 and certify pairs. Summing $\binom{\#\{p\in D:h_p>N-r\}}7$ over the interval gives exactly **4,631** minimum eight-prime lists with exterior prime above $N/2$, using **11** distinct exterior primes.

## Exhaustive enumeration

For each four-subset of $D$, the C++ scanner lists possible individual parts $a\le N/2$ by the Chinese remainder theorem. The allowed residues modulo $p$ are $0$ through $N\bmod p$. At each higher prime power $w$, a candidate with $a\bmod w>N\bmod w$ fails until the next multiple of $w$. This skips only forbidden candidates. Every row in increasing order has both first parts on the resulting list. For each admissible pair $a\le b\le N-a-b$ the test at every $w=p^i\le N$ is
$$
(a\bmod w)+(b\bmod w)\le N\bmod w.
$$
This is exactly the vanishing of that factorial-floor layer for the three parts. The independent Python check constructs the lowest residues from CRT idempotents, traverses each progression backwards, requires $a+b$ or its complement to occur among the small cuts, and tests actual base-$p$ additions. The two programs agree on all 1,716 four-subsets across the two weights. [C++ scanner](scan_endpoint_subset_rows.cpp) · [Python check](check_eight_prime_minima.py)

From the repository root, with a C++17 compiler and Python with `python-flint` installed, run

```sh
python3 check_six_prime_minimum.py
c++ -O3 -std=c++17 scan_endpoint_subset_rows.cpp -o /tmp/wasserman-four-subsets
python3 check_eight_prime_minima.py --scanner /tmp/wasserman-four-subsets --cross-all
```
