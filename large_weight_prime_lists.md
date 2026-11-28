# Direct common-prime lists at 7,364 selected weights

For each of the 7,364 weights from the [primitive-scale CRT search](primitive_scale_crt_search.md), the [table](large_weight_certificates.txt) lists primes such that **every three positive trinomial coefficients have a common prime on that line**. The largest weight is $10{,}940{,}257{,}758{,}976$. The certificates prove the assertion at their individual weights directly, using Kummer's criterion. They do not require the global reduction to the endpoint row.

At $N=4{,}217{,}316{,}279$, where a single coefficient has endpoint gcd below $N/8$, the list is
$$
          7,\quad 11,\quad 113,\quad 1321,\quad 1399,\quad
          1130951,\quad 4217316253 .
$$

| Listed primes | Weights |
|---:|---:|
| 1 | 10 |
| 4 | 7 |
| 5 | 765 |
| 6 | 1,882 |
| 7 | 3,396 |
| 8 | 1,073 |
| 9 | 195 |
| 10 | 33 |
| 11 | 3 |

The singleton certificates and the seven lists of size four are minimum. A universal prime must divide the endpoint $N(N-1)$ and it is universal precisely when $N$ has at most two base-$p$ digit units. Otherwise the digit units can be distributed among three nonempty parts without a carry. The generator checks every endpoint prime for universality. In the absence of a universal prime, one row per prime defeats any list of size at most three. Larger entries in the table give upper bounds.

## How the lists certify all triples

Let $r>N/2$ be the largest prime at most $N$. A row avoiding $r$ is exactly $(a,b,N-a-b)$ up to order, with $a,b>0$ and $a+b\le N-r$. If all three coefficients contain $r$, it is already a common divisor. Otherwise fix an avoiding row $R$. It suffices that every pair of arbitrary rows share a listed prime dividing $K_R$.

For most weights the remaining listed primes divide $N(N-1)$. Compute their supports on all possible $R$. Given such a support, distribute its primes arbitrarily between two putative rows that avoid them. The full powers of assigned weight primes have product $U$ and divide every part of their row. The normalized sum is $v=N/U$. At each assigned predecessor power $w$ the normalized parts have residues $v,0,0$ in some order. With $S$ the product of these powers, the CRT lists all triples of residues modulo $S$. Lift them to positive rows of the required sum and check all assigned bases. An allocation is impossible when at least one of its two sides has no such row. Exclude every allocation. This argument applies to subsets of endpoint primes even if some full powers from $N$ or $N-1$ are left unused.

Ten weights have single universal primes. At another 7,347 weights the preceding allocation test supplies the table directly. The generator conservatively skips individual side enumerations whose crude work estimate exceeds one million. Skipping can make it miss a shorter certificate. It accepts a list only after proving an exclusion for every required allocation.

## Seven weights using additional primes

At the seven remaining weights $N=3q$ or $4q$ for a large prime $q$. A row avoiding $q$ must respectively be $Q=(q,q,q)$ or $(q,q,2q)$. Its largest part is less than $r$, so a selection avoiding both $q$ and $r$ uses two distinct rows $Q,R$. Write $g=N-r$. If a full predecessor power $w\ge g$ has prime base $p$, every possible $R$ contains $p$ in its coefficient. Otherwise at least one of its two small parts would be a positive multiple of $w$, whereas their sum is at most $g$. Retain all such bases that also divide $K_Q$. The third row must avoid every retained base, or that base certifies the triple.

At six weights the product of the retained powers is $N-1$. At $542{,}587{,}407$ it is $(N-1)/2$. In each case the CRT and full carry test leave exactly **one** possible sorted third row. Factor the integers $N-j$ for $0\le j<g$. These contain every prime of every $K_R$ because its numerator is a falling factorial with at most $g$ terms. Test their primes on the unique third row, on $Q$, and on the exact coefficients of all $R$. Two additional primes cover all remaining choices in each case. The pairs are $\{2,79\}$ at $542{,}587{,}407$, $\{2,73\}$ at $2{,}222{,}751{,}452$, and $\{2,3\}$ at the other five weights. [Small-cofactor calculation](check_small_cofactor_certificates.py)

## Independent verification

The [independent verifier](verify_large_weight_certificates.py) checks primality and the supports on the $r$-avoiding rows by exact binomial products. For the two-row tests, it assigns CRT powers to the two smallest normalized parts via disjoint subsets and evaluates factorial floors. This enumerates the possibilities in a different order from the generator and uses exact valuations in place of its digit test. For the seven small-cofactor lists it directly imposes the supports common to the forced row and each near row, including any primes outside the endpoint. The table has one row per weight with primes in increasing order.

Run the [generator](check_large_weight_certificates.py) after the [primitive-scale search](primitive_scale_crt_search.md#reproduction), or run the standalone verification directly on the table. Both use Python and python-flint.

~~~sh
python3 verify_large_weight_certificates.py large_weight_certificates.txt
python3 check_large_weight_certificates.py .work/power_scale2000.rows .work/regenerated_lists.txt
cmp large_weight_certificates.txt .work/regenerated_lists.txt
~~~
